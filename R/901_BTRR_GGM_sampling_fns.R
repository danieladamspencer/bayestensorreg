#' Quadratic calculation
#'
#' @param beta_g region beta
#' @param W_g region omega
#'
#' @importFrom abind abind
#'
#' @return The results of the quadratic (D x Rank)
#' @keywords internal
BTRR_GGM_BtWB <- function(beta_g, W_g) {
  BtWB <- mapply(function(b, w) {
    bw_bind <- abind::abind(b, w, along = 3)
    apply(bw_bind, 2, function(bwb) {
      crossprod(bwb[, 1], diag(1 / bwb[, 2])) %*% bwb[, 1]
    })
  }, b = beta_g, w = W_g)
  BtWB <- as.matrix(BtWB)
  return(BtWB)
}

#' Check that the input data are compatible with the \code{BTRR_GGM} model function.
#'
#' @param input an input object to be checked
#'
#' @return If the input object passes all of the checks, then it is also returned.
#' @keywords internal
BTRR_GGM_check_data <- function(input) {
  is_list <- class(input) %in% c("list","TRR_GGM_data")
  if(!is_list) stop("The input data should be a list.")
  has_xY <- all(c("x","Y") %in% names(input))
  if(!has_xY) stop("The input should be a list with elements `x` and `Y`.")
  has_multiple_response_tensors <- length(input$Y) > 1
  if(!has_multiple_response_tensors) stop("There should be multiple response tensors within `input$Y`.")
  dim_Yg <- dim(input$Y[[1]])
  has_dim_DD <- length(dim_Yg) > 2
  if(!has_dim_DD) stop("Each response tensor should be an array of order at least 3 (at least one index for tensor location, an index for time, and an index for subject).")
  dims_equal <- all(sapply(input$Y,function(Yg) tail(dim(Yg),2) == tail(dim_Yg,2)))
  if(!dims_equal) stop("The last two dimensions (corresponding to time and subject) in the response tensors do not all match.")
  x_matches_Y <- all(dim(input$x) == tail(dim_Yg,2))
  if(!x_matches_Y) stop("The dimensions of the predictor do not match the final two dimensions of the response tensors.")
  class(input) <- "TRR_GGM_data"
  return(input)
}

#' Draw betas
#'
#' @param Y_g region response
#' @param x covariate
#' @param d_g region effect
#' @param beta_g activation effect
#' @param sig2y observation error
#' @param tau.g global variance parameter
#' @param phi.g local rank variance weighting
#' @param omega.g local variance
#' @param j dimension margin
#' @param r rank
#'
#' @return vector beta updates
#' @keywords internal
BTRR_GGM_draw_beta <- function(Y_g, x, d_g, beta_g, sig2y,
                               tau.g, phi.g, omega.g, j, r) {
  DD <- length(dim(Y_g)) - 2
  TT <- nrow(x)
  n <- ncol(x)
  if(length(phi.g) == 1){
    expected <- array(1,dim=head(dim(Y_g),-1)) %o% d_g
  }else{
    expected <-
      composeParafac(lapply(beta_g, function(beta_gj) {
        beta_gj[,-r, drop = FALSE]
      })) %o% x + array(1, dim = head(dim(Y_g), -1)) %o% d_g
  }
  y.til <- apply((Y_g - expected),c(j, (DD+1):(DD+2)), identity)
  y.til <- aperm(y.til, perm = c(2,1,3,4))
  B_g_noj <-
    composeParafac(sapply(beta_g[-j], function(beta_g_notj)
      beta_g_notj[, r, drop = F], simplify = F))
  var <- 1 / (sum(x ^ 2 %o% B_g_noj^ 2) / sig2y +
            (1 / (tau.g * phi.g[r]) / diag(omega.g[[j]][, r])))
  mean_beta <- var %*% apply(y.til,1,function(dim_margin){
    sum(sapply(seq(TT),function(each_time){
      sapply(seq(n),function(each_subject){
        x[each_time,each_subject] * B_g_noj * dim_margin[,each_time,each_subject]
      })
    })) / sig2y
  })
  beta_gjr <-
    rnorm(length(beta_g[[j]][,r]), mean_beta, sqrt(diag(var)))
  return(beta_gjr)
}

#' Draw d
#'
#' @param input the input data
#' @param params the list of parallel results
#' @param Sig The interregion precision
#' @param sig2y The observation variance
#'
#' @return the matrix d
#' @keywords internal
BTRR_GGM_draw_d <- function(input, params, Sig, sig2y) {
  DD <- length(dim(input$Y[[1]])) - 2
  TT <- nrow(input$x)
  p <- sapply(input$Y, function(Y_g) head(dim(Y_g),-2),simplify = F)
  y_hat <- mapply(function(y,parm){
    y - composeParafac(parm$betas) %o% input$x
  },y=input$Y,parm=params,SIMPLIFY = FALSE)

  inv_d_covar <- Sig + (TT*diag(sapply(p,prod)))^2 / sig2y
  d_covar <- chol2inv(chol(inv_d_covar))
  d_mean <- d_covar%*% t(sapply(y_hat, function(yh){apply(yh,DD+2,sum)}))

  d <- apply(d_mean,2,function(dm){
    dm + rnorm(length(dm)) %*% chol(d_covar)
  })
  return(d)
}

#' Update lambda
#'
#' @param a.lambda hyperparameter
#' @param b.lambda hyperparameter
#' @param beta_g region beta
#' @param phi.g region variance weights
#' @param tau.g region global variance
#'
#' @return update for lambda
#' @keywords internal
BTRR_GGM_draw_lambda <- function(a.lambda, b.lambda, beta_g, phi.g, tau.g) {
  Rank <- ncol(beta_g[[1]])
  sumabsb <- sapply(beta_g, function(beta_gj) {
    apply(beta_gj, 2, function(beta_gjr) {
      sum(abs(beta_gjr))
    })
  })
  if(Rank == 1){
    lambda.g <- sapply(sumabsb,function(sumabsb_j){
      rgamma(Rank,
             a.lambda + sapply(beta_g,nrow),
             b.lambda + (phi.g * tau.g) ^ (-.5) * sumabsb_j)
    })
    lambda.g <- t(lambda.g)
  }else{
    lambda.g <-
      apply(sumabsb, 2, function(sumabsb_j) {
        rgamma(Rank,
               a.lambda + sapply(beta_g,nrow),
               b.lambda + (phi.g * tau.g) ^ (-.5) * sumabsb_j)
      })
  }
  return(lambda.g)
}

#' Update local region variance parameters
#'
#' @param beta_g region location parameters
#' @param tau.g region global variance
#' @param phi.g region covariance weights
#' @param lambda.g region covariance shrinkage parameter
#'
#' @importFrom GIGrvg rgig
#'
#' @return an update for omega
#' @keywords internal
BTRR_GGM_draw_omega <- function(beta_g, tau.g, phi.g, lambda.g) {
  omega_g <- mapply(function(beta_gj, lambda_gj) {
    omega_gj <- mapply(function(beta_gjr, phi_gr) {
      chi <- beta_gjr ^ 2 / (tau.g * phi_gr)
      omega_gjr <- GIGrvg::rgig(nrow(beta_gj), 0.5, chi, lambda_gj)
    }, beta_gjr = split(beta_gj, col(beta_gj)), phi_gr = phi.g, SIMPLIFY = "array")
  }, beta_gj = beta_g, lambda_gj = split(lambda.g, row(lambda.g)),
  SIMPLIFY = FALSE)
  return(omega_g)
}

#' Draw Phi_g under the stick-breaking prior
#'
#' @param BtWB_g region normalized beta quadratic
#' @param Xi_g last draw region Xi
#' @param accept_g acceptance count for region g
#' @param cov_Metro_g metropolis covariance for region g
#' @param tau_g region tau
#' @param alpha region alpha
#'
#' @return List with a Phi update,  the acceptance number update, and the Xi
#'   update
#' @keywords internal
BTRR_GGM_draw_Phi <-
  function(BtWB_g,
           Xi_g,
           accept_g,
           cov_Metro_g,
           p_g,
           tau_g,
           alpha_g) {
    Rank = ncol(BtWB_g)
    old_Xi_g <- Xi_g
    if (Rank == 1) {
      phi.g <- 1
      accept <- 1
      old_Xi_g <- 1
    } else {
      accept <- accept_g
      new_Xi_g <-
        c(old_Xi_g + cov_Metro_g %*% rnorm(Rank - 1))
      while (length(new_Xi_g[new_Xi_g <= 0]) > 0) {
        new_Xi_g <- c(old_Xi_g + cov_Metro_g %*% rnorm(Rank - 1))
      }
      new_post_dens <- sum(sapply(seq(Rank - 1), function(cr) {
        stick_break_log_posterior(
          Xi_g = new_Xi_g,
          current_rank =  cr,
          BtWB_g =  BtWB_g,
          p_g = p_g,
          tau_g =  tau_g,
          alpha_g = alpha_g
        )
      }))
      old_post_dens <- sum(sapply(seq(Rank - 1), function(cr) {
        stick_break_log_posterior(
          Xi_g = old_Xi_g,
          current_rank =  cr,
          BtWB_g =  BtWB_g,
          p_g = p_g,
          tau_g =  tau_g,
          alpha_g = alpha_g
        )
      }))
      if (exp(new_post_dens - old_post_dens) > runif(1)) {
        old_Xi_g <- new_Xi_g
        accept <- accept + 1
      }
      phi.g <- stick_values(old_Xi_g)
    }
    return(list(
      phi.g = phi.g,
      accept.g = accept,
      Xi_g = old_Xi_g
    ))
  }

#' Draw Sigma and Upsilon
#'
#' @param d subject-region effects
#' @param g region
#' @param zeta scalar parameter
#' @param Sig Precision
#' @param Upsilon Precision parameter
#'
#' @importFrom statmod rinvgauss
#'
#' @return Updated Sig and Upsilon
#' @keywords internal
BTRR_GGM_draw_Sig_Upsilon <- function(d, g, zeta, Sig, Upsilon) {
  S <- tcrossprod(d)
  delta <- rgamma(1,ncol(d)/2 + 1,(S[g,g] + zeta)/2)
  S11inv <- chol2inv(chol(Sig[-g,-g]))
  varalpha <- chol2inv(chol((S[g,g] + zeta)*S11inv + diag(as.matrix(1/Upsilon[-g,g]))))
  meanalpha <- c(-varalpha%*%S[g,-g])
  s_alpha <- meanalpha + rnorm(nrow(d)-1) %*% chol(varalpha)
  Sig[g,-g] <- s_alpha
  Sig[-g,g] <- s_alpha
  Sig[g,g] <- delta + (s_alpha) %*% S11inv %*%t(s_alpha)
  if(g < nrow(d)){
    for(gg in (g+1):nrow(d)){
      Upsilon[gg,g] <- Upsilon[g,gg] <- 1/statmod::rinvgauss(1,sqrt(zeta^2 / Sig[g,gg]^2),zeta^2)
    }
  }
  return(list(Sig = Sig, Upsilon = Upsilon))
}

#' Update tau
#'
#' @param a.tau hyperparameter
#' @param b.tau hyperparameter
#' @param p_g region dimension
#' @param phi.g region rank variation weights
#' @param BtWB_g quadratic
#'
#' @importFrom GIGrvg rgig
#'
#' @return update for tau
#' @keywords internal
BTRR_GGM_draw_tau <- function(a.tau, b.tau, p_g, phi.g, BtWB_g) {
  Rank <- ncol(BtWB_g)
  if(Rank == 1){
    xi <- a.tau - sum(p_g) / 2
    chi <- sum(BtWB_g)
  }else{
    xi <- a.tau - ncol(BtWB_g) * sum(p_g) / 2
    chi <- sum(apply(BtWB_g, 1, sum) / phi.g)
  }
  tau.g <- GIGrvg::rgig(1, xi, chi, 2 * b.tau)
  return(tau.g)
}

#' Create a preallocated result list object
#'
#' @param p a list of length equal to the number of regions giving the response
#'   tensor dimensions for that region
#' @param Rank number of ranks in the model
#' @param n_iter number of iterations for the MCMC
#'
#' @return preallocated result list object
#' @keywords internal
BTRR_GGM_empty_results <- function(p, n, Rank, n_iter) {
  G <- length(p)
  check_D <- sapply(p,length)
  if(any(check_D !=  check_D[1]))
    stop("All of the regions should have the same dimension length.")
  DD <- unique(check_D)
  out <- list(
    betas = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(DD), function(j) {
          matrix(NA, p[[g]][j], Rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    W = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(DD), function(j) {
          matrix(NA, p[[g]][j], Rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    lambda = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(DD), function(j) {
          matrix(NA, DD, Rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    Phi = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        rep(NA,Rank)
      }, simplify = F)
    }, simplify = F),
    tau = matrix(NA, G, n_iter),
    d = array(NA, dim = c(G, n, n_iter)),
    zeta = numeric(n_iter),
    Sig = array(NA, dim = c(G, G, n_iter)),
    llik = vector(mode = "numeric",length = n_iter),
    sig2y = numeric(n_iter),
    alpha = numeric(n_iter),
    accept = vector(mode = "numeric",length = G)
  )
  return(out)
}

#' Extract the last iteration of the MCMC from a result file
#'
#' @param result_file A character string giving the relative file path to a
#'   desired result file.
#'
#' @importFrom statmod rinvgauss
#'
#' @return A list of the elements from the most recent MCMC draw
#' @keywords internal
BTRR_GGM_extract_last_iter <- function(result_file) {
  results_list <- readRDS(result_file)
  # Backward compatibility fix
  if("B" %in% names(results_list)) results_list$betas <- results_list$B
  out <- list()
  S <- length(results_list$betas)
  out$S <- S
  out$d <- results_list$d[,,S]
  out$betas <- results_list$betas[[S]]
  G <- length(out$betas)
  DD <- length(out$betas[[1]])
  Rank <- ncol(out$betas[[1]])
  if(G > 1){
    out$Sig <- results_list$Sig[,,S]
    out$zeta <- results_list$zeta[S]
    out$Upsilon <- matrix(1, G, G) ## Suggested by Wang
    diag(out$Upsilon) <- 0 ## Suggested by Wang
    for(g in seq(G)){
      if(g < G){
        for(gg in (g+1):G){
          out$Upsilon[gg,g] <- out$Upsilon[g,gg] <- 1/statmod::rinvgauss(1,sqrt(out$zeta^2 / out$Sig[g,gg]^2),out$zeta^2)
        }
      }
    }
  }
  out$tau <- results_list$tau[,S]
  # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
  out$alpha.g <-
    seq(Rank ^ (-DD), Rank ^ (-.1), length.out = 10)
  out$lambda <- results_list$lambda[[S]]
  out$W <- results_list$W[[S]]
  out$Phi <- results_list$Phi[[S]]
  out$sig2y <- results_list$sig2y[S]
  if(Rank > 1){
    out$Xi <- sapply(out$Phi,function(each_rank){
      out <- numeric(Rank - 1)
      out[1] <- each_rank[1]
      for(r in seq(2,Rank - 1)){
        out[r] <- each_rank[r] / (prod(head((1 - each_rank),r-1)))
      }
      return(out)
    })
  }else{
    out$Xi <- matrix(1,1,G)
  }
  out$cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(Rank - 1),simplify = FALSE)
  return(out)
}


#' Griddy Gibbs draw for alpha
#'
#' @param g region index
#' @param M number of samples per proposal
#' @param p region tensor dimensions
#' @param TT number of time steps
#' @param Rank model rank
#' @param cov_Metro Metropolis covariance
#' @param alpha.g grid of possible values for alpha
#' @param betas tensor decomposition components
#' @param W beta individual variances
#' @param Xi variance rank weight stick breaks
#' @param tau region-wide variance
#' @param a.tau hyperparameter
#' @param b.tau hyperparameter
#' @param d region-subject intercepts
#' @param input input data
#' @param sig2y observation variance
#'
#' @return a scalar value for alpha
#' @keywords internal
BTRR_GGM_griddy_Gibbs_alpha <- function(g, M, p, TT, Rank, cov_Metro, alpha.g, betas, W, Xi, tau, a.tau, b.tau, d, input, sig2y) {
  BtWB_g <- BTRR_GGM_BtWB(beta_g = betas[[g]],W_g = W[[g]])
  l.weights <- sapply(alpha.g, function(proposed) {
    bw <- mapply(function(b, w) {
      abind::abind(b, w, along = 3)
    },
    b = betas[[g]],
    w = W[[g]],
    SIMPLIFY = F)
    chi <- sapply(bw, function(each_dim) {
      apply(each_dim, 2, function(each_position) {
        t(each_position[, 1]) %*%
          diag(1 / each_position[, 2]) %*%
          each_position[, 1]
      })
    })
    chi <- apply(chi, 1, sum)
    ## Draw Phi proposals
    ##### Phi under a stick-breaking prior
    old_Xi_g <- Xi[,g]
    phi.l <- sapply(seq(M),function(m){
      new_Xi_g <- c(old_Xi_g + cov_Metro[[g]] %*%
                      rnorm(Rank - 1))
      while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
        new_Xi_g <- c(old_Xi_g + cov_Metro[[g]] %*%
                        rnorm(Rank - 1))
      }
      new_post_dens <- sum(sapply(seq(Rank - 1),function(cr){
        stick_break_log_posterior(Xi_g = new_Xi_g,
                                  current_rank = cr,
                                  BtWB_g =  BtWB_g,
                                  p_g = p[[g]],
                                  tau_g = tau[g],
                                  alpha_g = proposed)
      }))
      old_post_dens <- sum(sapply(seq(Rank - 1),function(cr){
        stick_break_log_posterior(Xi_g = old_Xi_g,
                                  current_rank = cr,
                                  BtWB_g =  BtWB_g,
                                  p_g = p[[g]],
                                  tau_g = tau[g],
                                  alpha_g = proposed)
      }))
      if(exp(new_post_dens - old_post_dens) > runif(1)) old_Xi_g <- new_Xi_g
      stick_values(old_Xi_g)
    })
    ## Draw tau proposals
    ### ANOTHER RANK 1 CHANGE
    if(Rank == 1){
      chi2 <- chi / phi.l
    }else{
      chi2 <- apply(phi.l, 2, function(each_proposal) {
        chi / each_proposal
      })
      chi2 <- colSums(chi2)
    }
    tau.l <- GIGrvg::rgig(M, a.tau - Rank * sum(p[[g]])/2, chi2, 2 * b.tau)
    refs <- list(phi = phi.l, tau = tau.l)
    ## Evaluate the densities
    lik.mean.tensor <-
      composeParafac(betas[[g]]) %o% input$x + array(1, dim = p[[g]]) %o% rep(1,TT) %o% d[g,]
    l.lik <-
      sum(dnorm(c(input$Y[[g]]), c(lik.mean.tensor), sqrt(sig2y), log = T)) # Log-likelihood
    l.bdens <-
      colSums(apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
        # Log prior density for all betas
        sapply(each_proposal[-1], function(each_rank_phi) {
          sum(unlist(sapply(bw, function(each_dim) {
            apply(each_dim, 2, function(each_rank_bw) {
              dnorm(
                each_rank_bw[, 1],
                0,
                each_proposal[1] * each_rank_phi * each_rank_bw[, 2],
                log = T
              )
            })
          })))
        })
      }))


    l.tau <-
      dgamma(refs$tau, a.tau, b.tau, log = T) # Log prior density for tau
    l.phi <-
      apply(refs$phi, 2, function(each_proposal) {
        lgamma(Rank * proposed) - Rank * lgamma(proposed) + sum((rep(proposed, Rank) - 1) * log(each_proposal))
      })
    # Log prior density for phi
    apply(cbind(l.phi, l.tau, l.bdens), 1, sum) + l.lik
  })
  mean.lweights <- apply(l.weights, 2, mean)
  weights <- exp(mean.lweights - max(mean.lweights))
  alpha <- sample(alpha.g, 1, prob = weights)
  return(alpha)
}

#' Get the log posterior under the stick-breaking prior for BTRR_GGM
#'
#' @param Xi_g region sticks
#' @param current_rank rank through the cycle
#' @param betas_g region beta
#' @param omega_g region omega
#' @param tau_g region tau
#' @param alpha_g region alpha
#'
#' @return The value for the log of the posterior density under the stick
#'   breaking prior
#' @keywords internal
stick_break_log_posterior <-
  function(Xi_g,
           current_rank,
           BtWB_g,
           p_g,
           tau_g,
           alpha_g) {
    model_rank <- length(Xi_g) + 1
    if (current_rank == model_rank)
      stop("This only needs to be done for r = 1,...,R-1!")
    part_a <- log(Xi_g[current_rank]) * (-sum(p_g) / 2)
    part_b <-
      (alpha_g - (model_rank - current_rank) * sum(p_g) / 2 - 1)*
      log(1 - Xi_g[current_rank])
    dim_sum_BWB <- apply(BtWB_g, 1, sum)
    part_c <- (1 / Xi_g[current_rank]) * dim_sum_BWB[current_rank]
    greater_ranks <-
      if (current_rank < model_rank - 1)
        seq(current_rank + 1, model_rank - 1)
    if (is.null(greater_ranks)) {
      part_d = 0
    } else{
      part_d <- sum(unlist(sapply(greater_ranks, function(each_rank) {
        (1 / (Xi_g[each_rank] * prod(1 - Xi_g[seq(each_rank - 1)]))) *
          dim_sum_BWB[each_rank]
      })))
    }
    part_e <-
      (1 / stick_values(Xi_g)[model_rank]) * dim_sum_BWB[model_rank]
    out <- part_a + part_b - (1 / tau_g) * (part_c + part_d + part_e)
    return(out)
  }
