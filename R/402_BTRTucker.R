#' Bayesian tensor regression with the Tucker decomposition
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param ranks A vector of length \code{length(dim(input$X)) - 1} giving the
#'   desired number of ranks per dimension of the tensor coefficient
#' @param n_iter (a scalar) the number of posterior samples desired
#' @param n_burn (a scalar) the number of posterior samples to discard as a
#'   burn-in
#' @param hyperparameters a list with the (scalar) elements \code{a.tau},
#'   \code{b.tau}, \code{a.lam}, \code{b.lam}, \code{nu}, \code{s_02},
#'   \code{a.sig}, \code{b.sig}, \code{Sig_0}, \code{mu_gam}, \code{alpha.grid},
#'    \code{a.u}, \code{b.u}, \code{a.z}, and/or \code{b.z},
#'   defining the values of the hyperparameters within the
#'   model. If \code{NULL}, then default values will be used. It is also
#'   possible to specify only a subset of the hyperparameters. The
#'   remaining hyperparameters are set to their default values.
#' @param save_dir (a character) A path to a directory in which the temporary
#'   results will be saved. Defaults to the current working directory. If
#'   \code{NULL}, no temporary saves are made.
#' @param CP Should the model be reduced to the CP decomposition? Default: FALSE
#'
#' @return A list with the posterior samples
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- BTRTucker(input)
#' }
BTRTucker <-
  function(input,
           ranks = rep(1,length(dim(input$X)) - 1),
           n_iter = 100,
           n_burn = 0,
           CP = FALSE,
           hyperparameters = NULL,
           save_dir = NULL) {
    # Logic checks for input
    if (n_burn > n_iter)
      stop("n_iter must be greater than n_burn.")
    if (any(ranks <= 0))
      stop("ranks must be a vector of positive integers.")
    if (length(input$eta) == 0)
      input$eta <- matrix(0, nrow = length(input$y), ncol = 1)
    # Pull input statistics
    Dim <- length(dim(input$X)) - 1
    if (length(ranks) != Dim)
      stop("The length of ranks should be equal to the dimension of each subject's tensor covariate.")
    if(CP & length(unique(ranks)) != 1)
      stop("If CP is TRUE, then all ranks must be equal.")
    n <- length(input$y)
    p <- head(dim(input$X), -1)

    # Set hyperparameters
    a.lam <- 3 # From Guhaniyogi et al. [2017]
    b.lam <- a.lam ^ (1 / (2 * Dim)) # From Guhaniyogi et al. [2017]
    nu <- 2 # From Guhaniyogi et al. [2017]
    s_02 <- -log(0.95)
    a.sig <- 3
    b.sig <- 20
    Sig_0 <-
      900 * diag(ncol(input$eta)) # From Guhaniyogi et al. [2017]
    mu_gam = rep(0, ncol(input$eta)) # From Guhaniyogi et al. [2017]
    alpha_grid <- sapply(ranks, function(r) {
      if (r == 1) {
        NULL
      } else{
        seq(r ^ (-Dim), r ^ (-.10), length.out = 10)
      }
    }, simplify = FALSE) # Based off of Guhaniyogi et al. [2017]
    a.tau <- 1 # This is what was in Shaan's code
    b.tau <- min(ranks) ^ (1 / Dim - 1) # From Guhaniyogi et al. [2017]
    a.u <- 3
    b.u <- a.u ^ (1 / (2 * Dim)) # Based off of Guhaniyogi et al. [2017]
    a.z <- 1
    b.z <-
      min(ranks) ^ (1 / Dim - 1) # Based off of Guhaniyogi et al. [2017]

    if(!is.null(hyperparameters))
      list2env(hyperparameters, envir = environment())

    # Preallocate storage
    results <-
      list(
        betas = mapply(
          function(ps, rr) {
            array(NA, dim = c(ps, rr, n_iter))
          },
          ps = p,
          rr = ranks,
          SIMPLIFY = FALSE
        ),
        W = mapply(
          function(ps, rr) {
            array(NA, dim = c(ps, rr, n_iter))
          },
          ps = p,
          rr = ranks,
          SIMPLIFY = FALSE
        ),
        tau = numeric(n_iter),
        lambda = sapply(ranks, function(r)
          matrix(NA, n_iter, r), simplify = FALSE),
        Pi = sapply(ranks, function(r)
          matrix(NA, n_iter, r), simplify = FALSE),
        gam = matrix(NA, n_iter, ncol(input$eta)),
        sig_y2 = numeric(n_iter),
        llik = numeric(n_iter),
        G = matrix(NA, prod(ranks), n_iter),
        z = numeric(n_iter),
        V = matrix(NA, prod(ranks), n_iter),
        U = matrix(NA, prod(ranks), n_iter)
      )

    # Set initial conditions
    V <- array(1, dim = ranks)
    U <- array(1, dim = ranks)
    W <- mapply(
      function(ps, r)
        matrix(1, ps, r),
      ps = p,
      r = ranks,
      SIMPLIFY = FALSE
    )
    lam <- sapply(ranks, rep, x = 1, simplify = FALSE)
    gam <- lm(input$y ~ -1 + input$eta)$coefficients
    gam <- ifelse(is.na(gam), 0, gam)
    y_til <- input$y - c(tcrossprod(t(gam), input$eta))
    avail_threads <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(avail_threads)
    B_init <- parallel::parApply(cl,input$X, seq(Dim), function(x, y_til) {
      out <- lm(y_til ~ -1 + x)$coefficients
      if(is.na(out)) out <- 0 # This is for some cases where there are values outside a mask
      return(out)
    }, y_til = y_til)
    parallel::stopCluster(cl)
    betas <- sapply(seq(Dim), function(d) {
      b_d <- svd(kFold(B_init,d), nu = ranks[d], nv = ranks[d])$u
      return(b_d)
    }, simplify = F)
    tau <- var(unlist(betas))
    vec_B_new <- Reduce(`%x%`,betas)
    vec_XB <- crossprod(vec_B_new,apply(input$X, Dim + 1, identity))
    if(!CP) {
      G_init <- lm(y_til ~ -1 + t(vec_XB))$coefficients
      G <- array(G_init, dim = ranks)
      z <- var(unlist(G))
    }
    if(CP) {
      ranks_along <- sapply(ranks, seq, simplify = FALSE)
      G <- Reduce(`%o%`, ranks_along)
      G <- as.numeric((G^(1/Dim)) %% 1 == 0)
      G <- array(G, dim = ranks)
    }
    sig_y2 <- var(y_til)
    # Begin MCMC
    start_MCMC <- proc.time()[3]
    for (s in seq(n_iter)) {
      tau <- BTRT_draw_tau(a.tau, b.tau, betas, W)
      lam <- BTRT_draw_lam(a.lam, b.lam, betas, tau)
      W <- BTRT_draw_omega(lam, betas, tau)
      y_til <- input$y - c(tcrossprod(t(gam), input$eta))
      for (j in seq(Dim)) {
        betas[[j]] <-
          BTRT_draw_Bj(
            y_til = y_til,
            X = input$X,
            tau = tau,
            betas = betas,
            W = W,
            sig_y2 = sig_y2,
            G = G,
            j = j
          )
      }

      if (!all(ranks == 1) | !CP) {
        z <- BTRT_draw_z(a.z, b.z, G, V)
        U <- BTRT_draw_U(a.u, b.u, G, z)
        V <- BTRT_draw_V(U, z, G)
        G <-
          BTRT_draw_G(
            y_til = y_til,
            betas = betas,
            X = input$X,
            z = z,
            V = V,
            sig_y2 = sig_y2
          )
      }
      # Draw sig_y2
      vecB <- c(Reduce(`%x%`, rev(betas)) %*% c(G))
      XB <- c(crossprod(matrix(c(input$X), ncol = tail(dim(input$X), 1)), vecB))
      y_til <- input$y -
        XB -
        c(input$eta %*% gam)
      sig_y2 <-
        BTRT_draw_sig_y2(a.sig = a.sig,
                        b.sig = b.sig,
                        y_til = y_til)
      # Draw gam
      if (!all(c(input$eta) == 0)) {
        y_til <- input$y - XB
        gam <- BTRT_draw_gam(input$eta, Sig_0, mu_gam, y_til, sig_y2)
      }
      # Find llik
      llik <-
        sum(dnorm(input$y,
                  XB  + c(tcrossprod(gam, input$eta)),
                  sqrt(sig_y2), log = TRUE))
      # Store the results
      results$betas[[s]] <- betas
      results$W[[s]] <- W
      results$tau[s] <- tau
      results$lambda[[s]] <- lam
      results$gam[s, ] <- gam
      results$sig_y2[s] <- sig_y2
      results$llik[s] <- llik
      results$G[,s] <- c(G)
      results$z[s] <- z
      results$V[,s] <- V
      results$U[,s] <- U
      # Status report
      if (s %% floor(n_iter / 10) == 0) {
        cat(
          "Iteration ",
          s,
          " of ",
          n_iter,
          ", llik = ",
          llik,
          ", time spent: ",
          proc.time()[3] - start_MCMC,
          ", est. time remaining:",
          ((proc.time()[3] - start_MCMC) / s) * (n_iter - s),
          "seconds \n"
        )
      }
      if (!is.null(save_dir) & s %% floor(n_iter / 10) == 0) {
        saveRDS(results, file.path(
          save_dir,
          paste0(
            "BTR_y_X_temp_rank_",
            paste0(ranks, collapse = ""),
            "_",
            format(Sys.Date(), "%Y%m%d"),
            ".rds"
          )
        ))
      }
    }
    results$total_time <- proc.time()[3] - start_MCMC
    # Remove burn-in
    if (n_burn > 0) {
      results$betas <- results$betas[-(1:n_burn)]
      results$W <- results$W[-(1:n_burn)]
      results$tau <- results$tau[-(1:n_burn)]
      results$lambda <- results$lambda[-(1:n_burn)]
      results$gam <- results$gam[-(1:n_burn), ]
      results$sig_y2 <- results$sig_y2[-(1:n_burn)]
      results$G <- results$G[,-seq(n_burn)]
      results$z <- results$z[-seq(n_burn)]
      results$V <- results$V[,-seq(n_burn)]
      results$U <- results$U[,-seq(n_burn)]
    }
    if(!all(ranks == 1)) results$G <- asplit(results$G, length(dim(results$G)))
    class(results) <- "BTRT_result"
    return(results)
  }
