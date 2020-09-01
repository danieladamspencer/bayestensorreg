#' Single-subject Bayesian Tensor Response Regression
#'
#' Performs the MCMC to draw from the posterior distribution for the model published by Guhaniyogi and Spencer (2020).
#'
#' @param input an object of class \code{TRR_data} or a list with elements \code{Y} (an array with dimensions \eqn{p_1\times \cdots \times p_D \times T \times 1}) and \code{x} (a \eqn{T \times 1} matrix).
#' @param n.iter (a scalar) the number of posterior samples desired
#' @param n.burn (a scalar) the number of posterior samples to discard as a burn-in
#' @param ranks (a positive integer) the number of ranks in the PARAFAC/CP tensor decomposition
#' @param hyperparameters a list with the (scalar) elements \code{a.tau}, \code{b.tau}, \code{a.lambda}, \code{b.lambda}, \code{a.epsilon}, and \code{b.epsilon} defining the values of the hyperparameters within the model. If \code{NULL}, then default values will be used.
#' @param save_after (an integer) An .rds file will be saved every \code{save_after} MCMC iterations with all of the results to that point (helpful for getting intermediate results in less time)
#' @param save_llik (a logical) Should the log-likelihood be calculated and saved at each iteration? Doing so comes at a cost, but can be used for model diagnostics. Defaults to \code{TRUE}.
#' @param save_dir (a character) A path to a directory in which the temporary results will be saved. Defaults to the current working directory.
#'
#' @import stats GIGrvg mvtnorm utils
#'
#' @return A list with the posterior samples
#' @export
#'
#' @examples
#' # Not run
#' # input <- TRR_simulated_data()
#' # results <- BTRR_single_subject(input)
BTRR_single_subject <-
  function(input,
           n.iter = 100,
           n.burn = 0,
           ranks = 1,
           hyperparameters = NULL,
           save_after = NULL,
           save_llik = TRUE,
           save_dir = ".") {
    #  > Gather metadata ----
    lowest_dim_margin <- min(head(dim(input$Y), -2))
    n <- tail(dim(input$Y), 1)
    p <- head(dim(input$Y), -2)
    D <- length(dim(input$Y)) - 2

  if(lowest_dim_margin < ranks) stop("The rank of your model cannot be larger than your smallest tensor dimension. Try a lower model rank.")
  if(is.vector(input$x)) input$x <- tcrossprod(input$x,rep(1,n))
  TT <- dim(input$x)[1] # Number of time steps

  # > Load necessary packages ----
  requireNamespace("GIGrvg") # Needed to draw lambda
  # > Hyperparameters ----
  a.epsilon <- a.lambda <- a.tau <- b.epsilon <- b.lambda <- b.tau <- NULL
  if(is.null(hyperparameters)) {
    hyperparameters <- list(
      a.tau = D - 1,
      b.tau = ranks^((1 / D) - 1), # These values are from Guhaniyogi et al [2017]
      a.lambda = 3,
      b.lambda = 3^(1/(2*D)), # These values are from Guhaniyogi et al [2017]
      a.epsilon = 1,
      b.epsilon = -log(0.95)
    )
  }

  list2env(hyperparameters,envir = environment())

  # > Set Storage ----
  results <-
    list(
      B = list(),
      W = list(),
      lambda = list(),
      Phi = list(),
      tau = vector(mode = "numeric",length = n.iter),
      llik = vector(mode = "numeric",length = n.iter),
      k = vector(mode = "numeric",length = n.iter),
      sigma_epsilon_sq = vector(mode = "numeric",length = n.iter),
      alpha = vector(mode = "numeric",length = n.iter)
    )

  # > Set Initials ----
  betas <-
    sapply(seq(D),function(each_dim){
      out <- matrix(rnorm(p[each_dim]*ranks),p[each_dim],ranks)
      return(out)
    },simplify = FALSE)
  tau <- 1
  alpha.g <-
    seq(ranks ^ (-2), ranks ^ (-.1), length.out = 10) # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
  lambda <- matrix(1, ranks, 2)
  W <-
    sapply(betas, function(each_dim) {
      array(1, dim = dim(each_dim))
    }, simplify = FALSE)
  Phi <- rep(1 / ranks, ranks)
  k <- 0
  sigma_epsilon_sq <- 1
  Sigma_AR <- toeplitz(k^seq(0,TT-1)) * sigma_epsilon_sq / (1 - k^2)
  if(ranks > 1){
    Xi <- rep(.6,ranks - 1)
    cov_Metro <- 0.001 * diag(ranks - 1)
    M <- 9
  }else{
    Xi <- 1
  }
  accept <- 0

  # > Run MCMC ----
  beginning_of_sampler <- proc.time()[3]
  for (s in 1:n.iter) {
    Cr <- mapply(function(b_j,W_j){
      mapply(function(b_jr,W_jr){
        crossprod(b_jr,diag(1/W_jr)) %*% b_jr
      },b_jr = split(b_j,col(b_j)), W_jr = split(W_j,col(W_j)))
    },b_j = betas, W_j = W) # R x D
    # >> Griddy-Gibbs ----
    if(ranks > 1){
      weights <- sapply(alpha.g,function(a){
        Xi_l <- sapply(seq(M),function(m){
          draw_Xi(Xi,cov_Metro,betas,W,tau,a)}) # R - 1 x M
        if(ranks == 2) {
          Phi_l <- sapply(Xi_l,stick_values)
        }else{
          Phi_l <- apply(Xi_l,2,stick_values)
        }
        # tau_l <- rgig(M, a.tau - ranks * sum(p)/2,
        #               apply(Phi_l,2,function(phi){
        #                 sum(rowSums(Cr)/phi)}),
        #               2*b.tau) # length M
        tau_l <- sapply(seq(M),function(m) {
          draw_tau(a.tau = a.tau,b.tau = b.tau,betas = betas,W = W,Phi = Phi_l)
        })
        if(ranks == 2) {
          ldXi <- sapply(Xi_l,function(xi) sum(dbeta(xi,1,a,log = TRUE)))
        }else{
          ldXi <- apply(Xi_l,2,function(xi) sum(dbeta(xi,1,a,log = TRUE)))
        }
        ldtau <- dgamma(tau_l,a*ranks,b.tau,log = T)
        tau_r <- apply(Phi_l,1,`*`,tau_l)
        ldbetas <- apply(tau_r,1,function(tr){sum(-(sum(p)/2)*log(tr) - rowSums(Cr)/tr)})
        out <- ldXi + ldtau + ldbetas
      })
      max_weight <- max(weights)
      alpha_probs <- apply(weights,2,function(w) mean(exp(w - max_weight)))
      alpha <- sample(alpha.g,size=1,prob = alpha_probs)
      Xi <- draw_Xi(Xi,cov_Metro,betas,W,tau,alpha)
      Phi <- stick_values(Xi)
    }else{
      alpha <- 1
      Phi <- 1
    }

    # >> Draw tau ----
    tau <- draw_tau(a.tau = a.tau,b.tau = b.tau,betas = betas,W = W,Phi = Phi)
    # >> Draw lambda ----
    sumabsb <- sapply(betas, function(beta_j) {
      apply(beta_j, 2, function(beta_jr) {
        sum(abs(beta_jr))
      })
    })
    if(ranks == 1){
      lambda <- sapply(sumabsb,function(each_dim){
        rgamma(ranks,
               a.lambda + p,
               b.lambda + (Phi * tau) ^ (-.5) * each_dim)
      })
      lambda <- t(lambda)
    }else{
      lambda <-
        apply(sumabsb, 2, function(each_dim) {
          rgamma(ranks,
                 a.lambda + p,
                 b.lambda + (Phi * tau) ^ (-.5) * each_dim)
        })
    }
    # >> Draw W ----
    W <- sapply(seq(D), function(each_dim) {
      sapply(seq(ranks), function(each_rank) {
        ch <-
          sapply(betas[[each_dim]][, each_rank], function(each_value) {
            each_value ^ 2 / (tau * Phi[each_rank])
          })
        GIGrvg::rgig(p[each_dim], 0.5, ch, lambda[each_rank, each_dim] ^ 2)
      }, simplify = T)
    }, simplify = F)

    # >> Draw betas ----
    for(r in seq(ranks)) {
      for(j in seq(D)) {
        betas[[j]][,r] <- draw_beta(Y = input$Y,x = input$x,betas,
                                    Sigma_AR = Sigma_AR,tau = tau,
                                    Phi = Phi,W = W,j = j,r = r)
      }
    }

    # Draw AR(1) parameters ----
    k <- draw_k(input$Y,input$x,betas,sigma_epsilon_sq)
    sigma_epsilon_sq <- draw_sigma_epsilon_sq(a.epsilon,b.epsilon,input$Y,input$x,betas,k)
    Sigma_AR <- toeplitz(k^seq(0,TT-1)) * sigma_epsilon_sq / (1 - k^2)

    # >> Get the log-likelihood ----
    if(save_llik == TRUE){
      requireNamespace("mvtnorm")
      Y_zero_mean <- input$Y - composeParafac(betas) %o% input$x
      llik <- sum(apply(Y_zero_mean,seq(D+2)[-(D+1)],function(Y_elli) {
        mvtnorm::dmvnorm(Y_elli,mean = rep(0,length(Y_elli)),sigma = Sigma_AR,log = TRUE)
      }))
      results$llik[s] <- llik
    }

    # >> Store the results ----
    results$B[[s]] <- betas
    results$W[[s]] <- W
    results$lambda[[s]] <- lambda
    results$Phi[[s]] <- Phi
    results$tau[s] <- tau
    results$alpha[s] <- alpha
    results$k[s] <- k
    results$sigma_epsilon_sq[s] <- sigma_epsilon_sq
    results$accept <- accept

    if(!is.null(save_after)){
      if(s %% save_after == 0){
        saveRDS(results, file = file.path(save_dir,paste0("400_",D,"D_rank_",ranks,"_first_",s,"_samples.rds")))
      }
    }

    if(s %% ceiling(n.iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Rank = ", ranks,
          " Iteration # ",s," of ",n.iter,
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - beginning_of_sampler, "  seconds #####\n",
          "##### Estimated time remaining: ", ((proc.time()[3] - beginning_of_sampler)/s)*(n.iter - s)," seconds #####\n"
        )
      )
    }
  } # End sampler

  if(n.burn > 0){
    results$B <- results$B[-(1:n.burn)]
    results$W <- results$W[-(1:n.burn)]
    results$lambda <- results$lambda[-(1:n.burn)]
    results$Phi <- results$Phi[-(1:n.burn)]
    results$tau <- results$tau[-(1:n.burn)]
  }
  results$accept <- results$accept / n.iter
  class(results) <- "BTRR_result"
  return(results)
}