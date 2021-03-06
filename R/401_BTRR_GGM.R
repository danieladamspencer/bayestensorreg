#' Bayesian Tensor Response Regression with Gaussian Graphical Model
#'
#' Performs the MCMC to draw from the posterior distribution for the model
#'   published by Spencer, Guhaniyogi, and Prado (2020).
#'
#' @param input an object of class \code{TRR_GGM_data} or a list with elements
#'   \code{Y} (a list of G arrays with dimensions
#'   \eqn{p_1\times \cdots \times p_D \times T \times n}) and
#'   \code{x} (a \eqn{T \times n} matrix).
#' @param n_iter (a scalar) the number of posterior samples desired
#' @param n_burn (a scalar) the number of posterior samples to discard as a
#'   burn-in
#' @param Rank (a positive integer) the rank for the PARAFAC/CP
#'   tensor decomposition
#' @param hyperparameters a list with named numbers containing at least one of
#'   the following: \code{a.tau}, \code{b.tau}, \code{a.lambda},
#'   \code{b.lambda}, \code{a.epsilon}, or \code{b.epsilon} defining the values
#'   of the hyperparameters within the model. If \code{NULL}, then default
#'   values will be used.
#' @param save_after (an integer) An .rds file will be saved every
#'   \code{save_after} MCMC iterations with all of the results to that point
#'   (helpful for getting intermediate results in less time)
#' @param save_llik (a logical) Should the log-likelihood be calculated and
#'   saved at each iteration? Doing so comes at a cost, but can be used for
#'   model diagnostics. Defaults to \code{TRUE}.
#' @param results_file (optional) The relative path to a result file. This is
#'   used to continue an MCMC chain on a set of data for which some results
#'   already exist.
#' @param num_cores The number of cores used for running the code in parallel
#' @param save_dir (a character) A path to a directory in which the temporary
#'   results will be saved. Defaults to the current working directory.
#'
#' @return A list object with the posterior draws from the MCMC chain.
#' @export
#'
#' @importFrom abind abind
#' @importFrom doParallel registerDoParallel
#' @import foreach
#'
#' @examples
#' \dontrun{
#' input <- TRR_GGM_simulated_data()
#' results <- BTRR_GGM(input)
#' }
BTRR_GGM <- function(input,
                     n_iter = 100,
                     n_burn = 0,
                     Rank = 1,
                     hyperparameters = NULL,
                     save_after = NULL,
                     save_llik = TRUE,
                     results_file = NULL,
                     num_cores = parallel::detectCores() - 2,
                     save_dir = ".") {
  input <- BTRR_GGM_check_data(input)
  if(class(input) != "TRR_GGM_data")
    stop("The input should be a list with class TRR_GGM_data")
  lowest_dim_margin <- min(sapply(input$Y, function(g)
    head(dim(g), -2)))
  if (Rank > lowest_dim_margin)
    stop(
      "The Rank should be less than or equal to the length\n
      of the smallest response tensor dimensiton."
    )
  if (length(input$Y) < 2)
    stop("If you only have one region of interest, use the BTRR_single_subject function.")
  dims_Y <- sapply(input$Y, dim)
  DD <- nrow(dims_Y) - 2
  if (any(c(dims_Y[(DD + 1), ], nrow(input$x)) != dims_Y[(DD + 1), 1]))
    stop("The length of time does not match up between the\n
         regions-of-interest and/or the design matrix.")
  if (any(dims_Y[(DD + 2), ] != dims_Y[(DD + 2), 1]))
    stop("The number of subjects does not match up \n
         between the regions-of-interest.")
  n <- tail(dim(input$Y[[1]]),1) # Number of subjects
  G <- length(input$Y) # Number of regions of interest
  p <- sapply(input$Y,
              function(Y_g) head(dim(Y_g),-2),
              simplify=FALSE)
  TT <- dim(input$x)[1] # Number of time steps

  # > Set up Parallelization ----
  requireNamespace("foreach",quietly = T)
  if(G > 1){
    num <- min(G, num_cores)
    if(num > 1){
      cl <- parallel::makeCluster(num)
      doParallel::registerDoParallel(cl)
      # parallel::clusterEvalQ(cl = cl, library(GIGrvg))
      # parallel::clusterEvalQ(cl = cl, library(doParallel))
      # parallel::clusterEvalQ(cl = cl, library(abind))
      # parallel::clusterEvalQ(cl = cl, library(truncnorm))
      # parallel::clusterEvalQ(cl = cl, library(bayestensorreg))
    }
  }

  # Hyperparameters
  a.tau = DD - 1
  b.tau = Rank^((1 / DD) - 1) # These values are from Guhaniyogi et al [2017]
  a.lambda = 3
  b.lambda = 3^(1/(2*DD)) # These values are from Guhaniyogi et al [2017]
  a.zeta = 1
  b.zeta = 0.01 # These values are from Wang [2012]
  a.sig = 1
  b.sig = -log(0.95) # These values are from Guhaniyogi et al [2017]

  if(!is.null(hyperparameters))
    list2env(hyperparameters, envir = environment())

  # Set Storage
  results <- BTRR_GGM_empty_results(p, n, Rank, n_iter)

  # > Set Initials ----
  # >> Continuing from other result files ----
  if(!is.null(results_file)){
    old_results <- BTRR_GGM_extract_last_iter(results_file)
    list2env(old_results, envir = environment())
  }else{
    d <- matrix(0,G,n)
    betas <- sapply(p, function(p_g) {
      sapply(p_g, function(p_gj) {
        matrix(rnorm(p_gj * Rank, sd = sd(unlist(input$Y))/100), p_gj, Rank)
      }, simplify = F)
    }, simplify = F)
    y_prime <- sapply(input$Y, function(Y_g) {
      lm_residuals <- apply(Y_g, seq(DD), function(y_gv) lm(y_gv ~ input$x)$residuals)
      lm_residuals <- aperm(lm_residuals, c(seq((DD+1))[-1],1))
      lm_residuals <- array(lm_residuals, dim = dim(Y_g))
      subject_sum_residuals <- apply(abs(lm_residuals), (DD+2), sum)
      return(subject_sum_residuals)
    }, simplify = T)
    Sig <- solve(crossprod(y_prime)/n)
    Sig_L <- diag(1,G*(G-1)/2)
    Upsilon <- matrix(1, G, G) ## Suggested by Wang
    diag(Upsilon) <- 0 ## Suggested by Wang
    zeta <- 3  ## Value used in Wang et al. [2012]
    tau <- rep(1, G) # Assuming unit value
    # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
    alpha.g <-
      seq(Rank ^ (-DD), Rank ^ (-.1), length.out = 10)
    lambda <-
      sapply(1:G, function(g) {
        matrix(1, Rank, DD)
      }, simplify = F)
    W <-
      lapply(betas, function(beta_g) {
        sapply(beta_g, function(each_dim) {
          array(1, dim = dim(each_dim))
        }, simplify = FALSE)
      }) # List of length 5, each element a list of length 2, each element a matrix dim p_j, Rank
    phi <-
      sapply(1:G, function(g) {
        rep(1 / Rank, Rank)
      }, simplify = F) # List of length G, each element a matrix of Rank rows by 2 columns
    sig2y <- 1 # Observational variance estimate
    if(Rank > 1){
      Xi <- sapply(seq(G),function(x) rep(.6,Rank - 1))
      if(!is.matrix(Xi)) Xi <- matrix(Xi,nrow = 1)
      cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(Rank - 1),simplify = FALSE)
    }else{
      Xi <- matrix(1,1,G)
    }
  }
  accept <- rep(0,G)

  # > Run MCMC ----
  beginning_of_sampler <- proc.time()[3]
  for (s in 1:n_iter) {
    # >> Griddy-Gibbs ----
    ## Almost everything is done separately based on the region of interest
    params <- foreach(g = seq(G),
                      .packages = c("GIGrvg","abind","bayestensorreg"),
                      .verbose = F, .errorhandling = "stop") %dopar% {
      # Set up to sample alpha IF RANK > 1
      if(Rank > 1 & s <= 100){
        alpha <- BTRR_GGM_griddy_Gibbs_alpha(g = g, M = 4, p = p, TT = TT,
                                             Rank = Rank, cov_Metro = cov_Metro,
                                             alpha.g = alpha.g, betas = betas,
                                             W = W, Xi = Xi, tau = tau,
                                             a.tau = a.tau, b.tau = b.tau,
                                             d = d, input = input, sig2y = sig2y)
      }else{
        alpha <- 0
      }
      BtWB_g <- BTRR_GGM_BtWB(betas[[g]], W[[g]])
      # >> Draw phi ----
      phi_draw <- BTRR_GGM_draw_Phi(BtWB_g = BtWB_g,
                                    Xi_g = Xi[,g],accept_g = results$accept[g],
                                    cov_Metro_g = cov_Metro[[g]],p_g = p[[g]],
                                    tau_g = tau[g], alpha_g = alpha)
      phi.g <- phi_draw$phi.g
      accept <- phi_draw$accept.g
      Xi.g <- phi_draw$Xi_g

      # >> Draw tau ----
      tau.g <-
        BTRR_GGM_draw_tau(
          a.tau = a.tau,
          b.tau = b.tau,
          p_g = p[[g]],
          phi.g = phi.g,
          BtWB_g = BtWB_g
        )
      # >> Draw lambda ----
      lambda.g <-
        BTRR_GGM_draw_lambda(
          a.lambda = a.lambda,
          b.lambda = b.lambda,
          beta_g = betas[[g]],
          phi.g = phi.g,
          tau.g = tau.g
        )
      # >> Draw omega ----
      omega.g <-
        BTRR_GGM_draw_omega(
          beta_g = betas[[g]],
          tau.g = tau.g,
          phi.g = phi.g,
          lambda.g = lambda.g
        )
      # >> Draw beta ----
      betas.g <- betas[[g]]
      # This next part has to be done sequentially, so the for loops are unavoidable
      for (each_dim in seq(DD)) {
        for (each_rank in 1:Rank) {
          betas.g[[each_dim]][,each_rank] <- BTRR_GGM_draw_beta(
            Y_g = input$Y[[g]],
            x = input$x,
            d_g = d[g, ],
            beta_g = betas.g,
            sig2y = sig2y,
            tau.g = tau.g,
            phi.g = phi.g,
            omega.g = omega.g,
            j = each_dim,
            r = each_rank
          )
          }
        }

      list(
        al = alpha,
        phi = phi.g,
        tau = tau.g,
        omega = omega.g,
        betas = betas.g,
        lambda = lambda.g,
        Xi = Xi.g,
        accept = accept
      )
    } # End dopar

    # >> Draw d ----
    d <- BTRR_GGM_draw_d(input, params, Sig, sig2y)
    # >> Draw Sig ----
    for(g in seq(G)) {
      Sig_Upsilon.g <- BTRR_GGM_draw_Sig_Upsilon(d,g,zeta,Sig,Upsilon)
      Sig <- Sig_Upsilon.g$Sig
      Upsilon <- Sig_Upsilon.g$Upsilon
    }
    if(s < 100) {
      zeta <- rgamma(1,a.zeta + G*(G+1)/2, b.zeta + sum(abs(Sig))/2)
    }

    # >> Draw sig2y ----
    if(G > 1){
      sum_sq_diff <- sum(unlist(mapply(function(y,parms,dd){
        (y - composeParafac(parms$betas) %o% input$x - array(1,dim=head(dim(y),-1))%o%dd)^2
      },y = input$Y,parms = params,dd = split(d,row(d)))))
    }else{
      sum_sq_diff <-
        sum(c(input$Y - composeParafac(params[[1]]$betas) %o% input$x)^2)
    }

    sig2y <- 1/rgamma(1,
                      a.sig + n*TT*sum(sapply(p,prod))/2,
                      b.sig + (1/2)*sum(sum_sq_diff))

    betas <- lapply(params, function(z) {
      z$betas
    })
    W <- lapply(params, function(z) {
      z$omega
    })
    lambda <- lapply(params, function(z) {
      z$lambda
    })
    phi <- lapply(params, function(z) {
      z$phi
    })
    tau <- sapply(params, function(z) {
      z$tau
    })

    if(Rank > 1){
      Xi <- sapply(params, function(z) z$Xi[drop = FALSE])
      if(!is.matrix(Xi)) Xi <- t(as.matrix(Xi))
    }
    accept_all <- sapply(params,function(z) z$accept)

    # >> Get the log-likelihood ----
    if(save_llik == TRUE){
      llik <- -0.5*log(2*pi*sig2y)*n*TT*G*sum(sapply(p,prod)) - 0.5/sig2y * sum(sum_sq_diff)
      results$llik[s] <- llik
    }


    results$betas[[s]] <- betas
    results$W[[s]] <- W
    results$lambda[[s]] <- lambda
    results$Phi[[s]] <- phi
    results$tau[, s] <- tau
    results$d[, , s] <- d
    results$Sig[, , s] <- Sig
    results$zeta[s] <- zeta
    results$sig2y[s] <- sig2y
    results$alpha[s] <- params[[1]]$al
    results$accept <- accept_all

    if(!is.null(save_after)){
      if(s %% save_after == 0){
        saveRDS(results,file = file.path(save_dir,paste0(DD,"D_rank_",Rank,"_first_",s,"_samples_",format(Sys.Date(),"%Y%m%d"),".rds")))
      }
    }

    if(s %% ceiling(n_iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Rank = ", Rank,
          " Iteration # ",s," of ",n_iter,
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - beginning_of_sampler, "  seconds #####\n",
          "##### Estimated time remaining: ", ((proc.time()[3] - beginning_of_sampler)/s)*(n_iter - s)," seconds #####\n"
        )
      )
    }

  } # End sampler

  if(n_burn > 0){
    results$betas <- results$betas[-(1:n_burn)]
    results$W <- results$W[-(1:n_burn)]
    results$lambda <- results$lambda[-(1:n_burn)]
    results$Phi <- results$Phi[-(1:n_burn)]
    results$tau <- results$tau[, -(1:n_burn)]
    results$d <- results$d[, , -(1:n_burn)]
    results$Sig <- results$Sig[, , -(1:n_burn)]
    results$zeta <- results$zeta[-(1:n_burn)]
    results$sig2y <- results$sig2y[-(1:n_burn)]
  }

  results$accept <- results$accept / n_iter
  results$total_time <- proc.time()[3] - beginning_of_sampler
  class(results) <- "BTRR_GGM_result"
  return(results)
} # End MCMC function
