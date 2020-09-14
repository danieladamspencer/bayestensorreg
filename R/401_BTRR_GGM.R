BTRR_GGM <- function(input,
                     n_iter = 100,
                     n_burn = 0,
                     rank = 1,
                     hyperparameters = NULL,
                     save_after = NULL,
                     save_llik = TRUE,
                     results_file = NULL,
                     save_dir = ".") {
  if(class(input) != "TRR_GGM_data")
    stop("The input should be a list with class TRR_GGM_data")
  lowest_dim_margin <- min(sapply(input$Y, function(g)
    head(dim(g), -2)))
  if (rank > lowest_dim_margin)
    stop(
      "The rank should be less than or equal to the length\n
      of the smallest response tensor dimensiton."
    )
  if (length(input$Y) < 2)
    stop("If you only have one region of interest, use the BTRR_single_subject function.")
  dims_Y <- sapply(input$Y, dim)
  D <- nrow(dims_Y) - 2
  if (any(c(dims_Y[(D + 1), ], nrow(input$x)) != dims_Y[(D + 1), 1]))
    stop("The length of time does not match up between the\n
         regions-of-interest and/or the design matrix.")
  if (any(dims_Y[(D + 2), ] != dims_Y[(D + 2), 1]))
    stop("The number of subjects does not match up \n
         between the regions-of-interest."
    )
  n <- tail(dim(input$Y[[1]]),1) # Number of subjects
  G <- length(input$Y) # Number of regions of interest
  p <- sapply(input$Y,
              function(Y_g) head(dim(Y_g),-2),
              simplify=FALSE)
  TT <- dim(input$x)[1] # Number of time steps

  # > Load necessary packages ----
  # require(GIGrvg)
  # require(mvnfast)
  # require(doParallel)
  # require(abind)
  # require(dlm)
  # require(statmod)
  # require(truncnorm)
  # require(Rcpp)
  # require(RcppArmadillo)

  # > Functions ----
  # source("R/900_misc.R")

  # > Set up Parallelization ----
  # if(G > 1){
  #   num <- min(G, detectCores() - 1)
  #   if(num > 1){
  #     cl <- makeCluster(num)
  #     registerDoParallel(cl)
  #     clusterEvalQ(cl = cl, library(GIGrvg))
  #     clusterEvalQ(cl = cl, library(doParallel))
  #     clusterEvalQ(cl = cl, library(abind))
  #     clusterEvalQ(cl = cl, library(truncnorm))
  #     clusterEvalQ(cl = cl, library(Rcpp))
  #     clusterEvalQ(cl = cl, library(RcppArmadillo))
  #     clusterEvalQ(cl = cl, source("R/900_misc.R"))
  #   }
  # }

  # Hyperparameters
  a.tau = D - 1
  b.tau = rank^((1 / D) - 1) # These values are from Guhaniyogi et al [2017]
  a.lambda = 3
  b.lambda = 3^(1/(2*D)) # These values are from Guhaniyogi et al [2017]
  a.zeta = 1
  b.zeta = 0.01 # These values are from Wang [2012]
  a.sig = 1
  b.sig = -log(0.95) # These values are from Guhaniyogi et al [2017]

  if(!is.null(hyperparameters))
    list2env(hyperparameters, envir = environment())

  # Set Storage
  results <- BTRR_GGM_empty_results(p, rank, n_iter)

  # > Set Initials ----
  # >> Continuing from other result files ----
  if(!is.null(results_file)){
    old_results <- BTRR_GGM_extract_last_iter(results_file)
    list2env(old_results, envir = environment())
  }else{
    # >> Using best-guess initial values ----
    d <- matrix(0,G,n)
    betas <- sapply(p, function(p_g) {
      sapply(p_g, function(p_gj) {
        matrix(rnorm(p_gj * rank, sd = sd(unlist(input$Y))/100), p_gj, rank)
      }, simplify = F)
    }, simplify = F)
    y_prime <- sapply(input$Y, function(Y_g) {
      lm_residuals <- apply(Y_g, seq(D), function(y_gv) lm(y_gv ~ input$x)$residuals)
      lm_residuals <- aperm(lm_residuals, c(2,3,1))
      lm_residuals <- array(lm_residuals, dim = dim(Y_g))
      subject_sum_residuals <- apply(abs(lm_residuals), (D+2), sum)
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
      seq(rank ^ (-D), rank ^ (-.1), length.out = 10)
    lambda <-
      sapply(1:G, function(g) {
        matrix(1, rank, 2)
      }, simplify = F)
    W <-
      lapply(betas, function(beta_g) {
        sapply(beta_g, function(each_dim) {
          array(1, dim = dim(each_dim))
        }, simplify = FALSE)
      }) # List of length 5, each element a list of length 2, each element a matrix dim p_j, rank
    phi <-
      sapply(1:G, function(g) {
        rep(1 / rank, rank)
      }, simplify = F) # List of length G, each element a matrix of rank rows by 2 columns
    sig2y <- 1 # Observational variance estimate
    if(rank > 1){
      Xi <- sapply(seq(G),function(x) rep(.6,rank - 1))
      if(!is.matrix(Xi)) Xi <- matrix(Xi,nrow = 1)
      cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(rank - 1),simplify = FALSE)
    }else{
      Xi <- matrix(1,1,G)
    }
  }
  accept <- rep(0,G)

  # > Run MCMC ----
  beginning_of_sampler <- proc.time()[3]
  # pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
  for (s in 1:n_iter) {
    # >> Griddy-Gibbs ----
    ## Almost everything is done separately based on the region of interest
    params <- foreach(g = seq(G)) %dopar% {
      # Set up to sample alpha IF RANK > 1
      if(rank > 1 & s <= 100){
        M = 4 # Number of reference sets per grid value of alpha
        l.weights <- sapply(alpha.g, function(proposed) {
          bw <- mapply(function(b, w) {
            abind(b, w, along = 3)
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
          ### INCLUDE RANK 1 change:
          if(rank == 1){
            chi <- sum(chi)
          }else{
            chi <- apply(chi, 1, sum)
          }
          # This helps to debug when values go to infinity. Eventually, a workaround
          # should make this obsolete
          for (abc in 1:length(chi)) {
            if (is.infinite(chi[abc]))
              browser()
            if (is.nan(chi[abc]))
              browser()
          }
          ## Draw Phi proposals

          ##### Phi under a stick-breaking prior
          old_Xi_g <- Xi[,g]
          phi.l <- sapply(seq(M),function(m){
            new_Xi_g <- c(old_Xi_g + cov_Metro[[g]] %*%
                            rnorm(rank - 1))
            while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
              new_Xi_g <- c(old_Xi_g + cov_Metro[[g]] %*%
                              rnorm(rank - 1))
            }
            new_post_dens <- sum(sapply(seq(rank - 1),function(cr){
              stick_break_log_posterior(new_Xi_g,cr, betas[[g]],W[[g]],tau[g],proposed)
            }))
            old_post_dens <- sum(sapply(seq(rank - 1),function(cr){
              stick_break_log_posterior(old_Xi_g,cr, betas[[g]],W[[g]],tau[g],proposed)
            }))
            if(exp(new_post_dens - old_post_dens) > runif(1)) old_Xi_g <- new_Xi_g
            stick_values(old_Xi_g)
          })
          ## Draw tau proposals
          ### ANOTHER RANK 1 CHANGE
          if(rank == 1){
            chi2 <- chi / phi.l
          }else{
            chi2 <- apply(phi.l, 2, function(each_proposal) {
              chi / each_proposal
            })
            chi2 <- colSums(chi2)
          }
          tau.l <- rgig(M, a.tau - rank * sum(p[[g]])/2, chi2, 2 * b.tau)
          # tau.l <- rep(1,M)
          refs <- list(phi = phi.l, tau = tau.l)
          ## Evaluate the densities
          lik.mean.tensor <-
            composeParafac(betas[[g]]) %o% input$x + array(1, dim = p[[g]]) %o% rep(1,TT) %o% d[g,]
          # l.lik <-
          #   sum(dnorm(unlist(yg[[g]]), c(lik.mean.tensor), sqrt(sig2y), log = T)) # Log-likelihood
          l.lik <-
            sum(dnorm(c(input$Y[[g]]), c(lik.mean.tensor), sqrt(sig2y), log = T)) # Log-likelihood
          if(rank == 1){
            l.bdens <- apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
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
            })
          }else{
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
          }

          l.tau <-
            dgamma(refs$tau, a.tau, b.tau, log = T) # Log prior density for tau

          ### RANK 1 CHANGE:
          if(rank == 1){
            l.phi <- sapply(refs$phi,function(each_proposal){
              lgamma(rank * proposed) - rank * lgamma(proposed) + sum((rep(proposed, rank) - 1) * log(each_proposal))
            })
          }else{
            l.phi <-
              apply(refs$phi, 2, function(each_proposal) {
                lgamma(rank * proposed) - rank * lgamma(proposed) + sum((rep(proposed, rank) - 1) * log(each_proposal))
              })
          }
          # Log prior density for phi
          apply(cbind(l.phi, l.tau, l.bdens), 1, sum) + l.lik
        })
        mean.lweights <- apply(l.weights, 2, mean)
        weights <- exp(mean.lweights - max(mean.lweights))
        alpha <- sample(alpha.g, 1, prob = weights)
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
      Xi[,g] <- phi_draw$Xi_g

      # bg <- betas[[g]]
      # wg <- W[[g]]
      # ch <- mapply(function(b, w) {
      #   apply(abind(b, w, along = 3), 2, function(each_rank) {
      #     crossprod(each_rank[, 1], diag(1 / each_rank[, 2])) %*% each_rank[, 1]
      #   })
      # }, b = bg, w = wg)
      #
      # ##### Phi under a stick-breaking prior
      # old_Xi_g <- Xi[,g]
      # if(rank == 1){
      #   phi.g <- 1
      # }else{
      #   accept <- results$accept[g]
      #   new_Xi_g <- c(old_Xi_g + cov_Metro[[g]]%*%rnorm(rank - 1))
      #   while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
      #     new_Xi_g <- c(old_Xi_g + cov_Metro[[g]]%*%rnorm(rank - 1))
      #   }
      #   new_post_dens <- sum(sapply(seq(rank - 1),function(cr){
      #     stick_break_log_posterior(new_Xi_g,cr, betas[[g]],W[[g]],tau[g],alpha)
      #   }))
      #   old_post_dens <- sum(sapply(seq(rank - 1),function(cr){
      #     stick_break_log_posterior(old_Xi_g,cr, betas[[g]],W[[g]],tau[g],alpha)
      #   }))
      #   if(exp(new_post_dens - old_post_dens) > runif(1)){
      #     old_Xi_g <- new_Xi_g
      #     accept <- accept + 1
      #   }
      #   phi.g <- stick_values(old_Xi_g)
      # }

      # >> Draw tau ----
      tau.g <-
        BTRR_GGM_draw_tau(
          a.tau = a.tau,
          b.tau = b.tau,
          p_g = p[[g]],
          phi.g = phi.g,
          BtWB_g = BtWB_g
        )
      # xi <- a.tau - rank * sum(p[[g]])/2
      # ### RANK 1 ADJUSTMENT
      # if(rank == 1){
      #   chi <- sum(ch)
      # }else{
      #   chi <- sum(apply(ch, 1, sum) / phi.g)
      # }
      # tau.g <- rgig(1, xi, chi, 2 * b.tau)

      # >> Draw lambda ----
      lambda.g <-
        BTRR_GGM_draw_lambda(
          a.lambda = a.lambda,
          b.lambda = b.lambda,
          beta_g = betas[[g]],
          phi.g = phi.g,
          tau.g = tau.g
        )
      # sumabsb <- sapply(bg, function(each_dim) {
      #   apply(each_dim, 2, function(each_rank) {
      #     sum(abs(each_rank))
      #   })
      # })
      # if(rank == 1){
      #   lambda.g <- sapply(sumabsb,function(each_dim){
      #     rgamma(rank,
      #            a.lambda + p[[g]],
      #            b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
      #   })
      #   lambda.g <- t(lambda.g)
      # }else{
      #   lambda.g <-
      #     apply(sumabsb, 2, function(each_dim) {
      #       rgamma(rank,
      #              a.lambda + p[[g]],
      #              b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
      #     })
      # }

      # >> Draw omega ----
      omega.g <-
        BTRR_GGM_draw_omega(
          beta_g = betas[[g]],
          tau.g = tau.g,
          phi.g = phi.g,
          lambda.g = lambda.g
        )
      # omega.g <- sapply(seq(D), function(each_dim) {
      #   sapply(seq(rank), function(each_rank) {
      #     ch <-
      #       sapply(bg[[each_dim]][, each_rank], function(each_value) {
      #         each_value ^ 2 / (tau.g * phi.g[each_rank])
      #       })
      #     rgig(p[[g]][each_dim], 0.5, ch, lambda.g[each_rank, each_dim])
      #   }, simplify = T)
      # }, simplify = F)

      # >> Draw beta ----
      # This next part has to be done sequentially, so the for loops are unavoidable
      for (each_dim in seq(D)) {
        for (each_rank in 1:rank) {
          if(rank == 1){
            expected <- array(1,dim=c(p[[g]],TT)) %o% d[g,]
          }else{
            expected <- composeParafac(lapply(bg, function(each_dim_I) {each_dim_I[, -each_rank,drop = FALSE]})) %o% input$x + array(1,dim=c(p[[g]],TT)) %o% d[g,]
          }
          if(G > 1){
            y.til <- (input$Y[[g]] - expected) %>%
              apply((D+1):(D+2),mode_k_matriz,each_dim) %>%
              array(dim=c(p[[g]][each_dim],prod(p[[g]][-each_dim]),TT,n))
          }else{
            y.til <- (input$Y - expected) %>%
              apply((D+1):(D+2),mode_k_matriz,each_dim) %>%
              array(dim=c(p[[g]][each_dim],prod(p[[g]][-each_dim]),TT,n))
          }
          betas_by_rank <-
            sapply(seq(rank), function(rr) {
              sapply(seq(D), function(dd) {
                bg[[dd]][, rr]
              },simplify = FALSE)
            },simplify = FALSE) # This restructuring allows for the calculation
          vec_outer_other_betas <- c(Reduce(outer,betas_by_rank[[each_rank]][-each_dim]))
          var <-  1 /  (sum(input$x^2 %o% vec_outer_other_betas^2) / sig2y + (1 / (tau.g * phi.g[each_rank]) / diag(omega.g[[each_dim]][,each_rank]) ))
          mean_beta <- var %*% apply(y.til,1,function(dim_margin){
            sum(sapply(seq(TT),function(each_time){
              sapply(seq(n),function(each_subject){
                input$x[each_time,each_subject] * vec_outer_other_betas * dim_margin[,each_time,each_subject]
              })
            })) / sig2y
          })
          bg[[each_dim]][, each_rank] <-
            rnorm(p[[g]][each_dim], mean_beta, sqrt(diag(var)))
          if (each_dim > 1 && each_rank > 1) {
            bg[[each_dim]][1, each_rank] <-
              rtruncnorm(
                1,
                b = bg[[each_dim]][1, (each_rank - 1)],
                mean = (mean_beta)[1],
                sd = sqrt(var)[1]
              )
          }
        }
      }

      list(
        al = alpha,
        phi = phi.g,
        tau = tau.g,
        omega = omega.g,
        betas = bg,
        lambda = lambda.g,
        Xi = old_Xi_g[drop = FALSE],
        accept = accept
      )
    } # End dopar

    # >> Draw d ----
    y_hat <- mapply(function(y,parm){
      y - composeParafac(parm$betas) %o% input$x
    },y=input$Y,parm=params,SIMPLIFY = FALSE)

    inv_d_covar <- Sig + (TT*diag(sapply(p,prod)))^2 / sig2y
    #  #####################################################
    #  # Inverse Wishart prior on Sigma
    # (Or when the prior is applied to the covariance rather than the precision)
    # inv_d_covar <- chol2inv(chol(Sig)) + (TT*diag(sapply(p,prod)))^2 / sig2y
    #  #####################################################
    d_covar <- chol2inv(chol(inv_d_covar))
    d_mean <- d_covar%*% t(sapply(y_hat, function(yh){apply(yh,D+2,sum)}))

    d <- apply(d_mean,2,function(dm){
      dm + rnorm(length(dm)) %*% chol(d_covar)
    })

    # >> Draw Sig ----
    matd <- t(d)
    S <- crossprod(matd,matd)
    # S <- chol2inv(chol(S)) # Will this work if I'm looking at a covariance rather than a precision?
    for(g in 1:G){
      delta <- rgamma(1,n/2 + 1,(S[g,g] + zeta)/2)
      S11inv <- chol2inv(chol(Sig[-g,-g]))
      varalpha <- chol2inv(chol((S[g,g] + zeta)*S11inv + diag(as.matrix(1/Upsilon[-g,g]))))
      meanalpha <- c(-varalpha%*%S[g,-g])
      s_alpha <- meanalpha + rnorm(G-1) %*% chol(varalpha)
      Sig[g,-g] <- s_alpha
      Sig[-g,g] <- s_alpha
      Sig[g,g] <- delta + (s_alpha) %*% S11inv %*%t(s_alpha)
      if(g < G){
        for(gg in (g+1):G){
          Upsilon[gg,g] <- Upsilon[g,gg] <- 1/rinvgauss(1,sqrt(zeta^2 / Sig[g,gg]^2),zeta^2)
        }
      }
    }
    if(s < 100) {
      zeta <- rgamma(1,a.zeta + G*(G+1)/2, b.zeta + sum(abs(Sig))/2)
    }

    ################################################
    # Elliptical Slice Sampler for Sigma
    # nu <- rnorm(sum(seq(G)))
    # u <- runif(1)
    # upper_tri_Sig <- c(Sig[upper.tri(Sig,diag = TRUE)])
    # log_L_Sig <- sum(sapply(seq(n),function(i){
    #   mvnfast::dmvn(t(d[,i]),
    #                 d_mean[,i],
    #                 Sig,
    #                 log = TRUE,isChol = FALSE)
    # })) +
    #   mvnfast::dmvn(t(Sig[upper.tri(Sig,diag = TRUE)]),
    #                 rep(0,sum(seq(G))),
    #                 diag(1,sum(seq(G))),
    #                 log = TRUE,isChol = TRUE)
    # log_L_threshold <- log_L_Sig + log(u)
    # theta <- runif(1,0,2*pi)
    # theta_min <- theta - 2*pi
    # theta_max <- theta
    # log_L_Sig_prime <- -Inf
    # cat("+++++ Starting Iteration +++++")
    # sig_step <- 1
    # while(log_L_Sig_prime <= log_L_threshold){
    #   Sig_prime <- upper_tri_Sig*cos(theta) + nu*sin(theta)
    #   chol_Sig_prime <- matrix(0,G,G)
    #   chol_Sig_prime[upper.tri(chol_Sig_prime,diag = TRUE)] <- Sig_prime
    #   log_L_Sig_prime <- sum(sapply(seq(n),function(i){
    #     mvnfast::dmvn(t(d[,i]),
    #                   d_mean[,i],
    #                   tcrossprod(chol_Sig_prime),
    #                   log = TRUE,isChol = FALSE)
    #   })) +
    #     mvnfast::dmvn(t(Sig_prime),
    #                   rep(0,sum(seq(G))),
    #                   diag(1,sum(seq(G))),
    #                   log = TRUE,isChol = TRUE)
    #   if(theta < 0){
    #     theta_min <- theta
    #   }else{
    #     theta_max <- theta
    #   }
    #   theta <- runif(1,theta_min,theta_max)
    #   cat("----- Ended iteration",sig_step,":",theta,"in [",theta_min,",",theta_max,"] ----- \n")
    #   sig_step <- sig_step + 1
    #   if(sig_step == 50) return(results)
    #   stopifnot(sig_step < 50)
    # }
    # Sig <- tcrossprod(chol_Sig_prime)
    ################################################

    ################################################
    # Inverse Wishart Prior on Sigma
    # inv_Sig <- MCMCpack::rwish(n+nu,VV + tcrossprod(d))
    # Sig <- chol2inv(chol((inv_Sig + t(inv_Sig))/2 + 1e-8))
    ################################################


    betas <- lapply(params, function(z) {
      z$betas
    })

    # >> Draw sig2y ----
    if(G > 1){
      sum_sq_diff <- mapply(function(y,parms,dd){
        (y - composeParafac(parms$betas) %o% input$x - array(1,dim=head(dim(y),-1))%o%dd)^2
      },y = input$Y,parms = params,dd = split(d,row(d))) %>% unlist %>% sum
    }else{
      sum_sq_diff <-
        (input$Y - composeParafac(params[[1]]$betas) %o% input$x)^2 %>%
        c %>%
        sum
    }

    sig2y <- 1/rgamma(1,
                      a.sig + n*TT*sum(sapply(p,prod))/2,
                      b.sig + (1/2)*sum(sum_sq_diff))

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

    if(rank > 1){
      Xi <- sapply(params, function(z) z$Xi[drop = FALSE])
      if(!is.matrix(Xi)) Xi <- t(as.matrix(Xi))
    }
    accept_all <- sapply(params,function(z) z$accept)

    # >> Get the log-likelihood ----
    if(save_llik == TRUE){
      llik <- -0.5*log(2*pi*sig2y)*n*TT*G*sum(sapply(p,prod)) - 0.5/sig2y * sum(sum_sq_diff)
      results$llik[s] <- llik
    }


    results$B[[s]] <- betas
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

    # save(results,s,file="C:/Users/Dan/Desktop/temp.Rdata")

    # if(save_results & s %% 10 == 0) save(results,file = paste0("intermediate_results",format(Sys.Date(),"%Y%m%d"),".RData"))

    if(!is.null(save_after)){
      if(s %% save_after == 0){
        saveRDS(results,file = file.path(save_dir,paste0("405_",D,"D_rank_",rank,"_first_",s,"_samples_",format(Sys.Date(),"%Y%m%d"),".rds")))
      }
    }

    # setTxtProgressBar(pb, s)
    if(s %% ceiling(n_iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Rank = ", rank,
          " Iteration # ",s," of ",n_iter,
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - beginning_of_sampler, "  seconds #####\n",
          "##### Estimated time remaining: ", ((proc.time()[3] - beginning_of_sampler)/s)*(n_iter - s)," seconds #####\n"
        )
      )
    }

  } # End sampler

  if(n_burn > 0){
    results$B <- results$B[-(1:n_burn)]
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
  return(results)
} # End MCMC function
