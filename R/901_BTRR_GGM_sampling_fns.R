#' Quadratic calculation
#'
#' @param beta_g region beta
#' @param W_g region omega
#'
#' @return The results of the quadratic (D x rank)
#' @keywords internal
BTRR_GGM_BtWB <- function(beta_g, W_g) {
  BtWB <- mapply(function(b, w) {
    bw_bind <- abind::abind(b, w, along = 3)
    apply(bw_bind, 2, function(bwb) {
      crossprod(bwb[, 1], diag(1 / bwb[, 2])) %*% bwb[, 1]
    })
  }, b = beta_g, w = W_g)
  return(BtWB)
}

BTRR_GGM_draw_beta <- function(Y_g, x, d_g, beta_g, sig2y,
                               tau.g, phi.g, omega.g, j, r) {
  if(length(phi.g) == 1){
    expected <- array(1,dim=head(dim(Y_g),-1)) %o% d_g
  }else{
    expected <-
      composeParafac(lapply(beta_g, function(beta_gj) {
        beta_gj[,-r, drop = FALSE]
      })) %o% x + array(1, dim = head(dim(Y_g), -1)) %o% d_g
  }
  y.til <- apply((Y_g - expected),c(j, (D+1):(D+2)), identity) #%>%
  y.til <- aperm(y.til, perm = c(2,1,3,4))
  B_g_noj <-
    composeParafac(sapply(beta_g[-j], function(beta_g_notj)
      beta_g_notj[, r, drop = F], simplify = F))
  var <- 1 / (sum(x ^ 2 %o% vec_outer_other_betas ^ 2) / sig2y +
            (1 / (tau.g * phi.g[r]) / diag(omega.g[[j]][, r])))
  mean_beta <- var %*% apply(y.til,1,function(dim_margin){
    sum(sapply(seq(TT),function(each_time){
      sapply(seq(n),function(each_subject){
        input$x[each_time,each_subject] * vec_outer_other_betas * dim_margin[,each_time,each_subject]
      })
    })) / sig2y
  })
  beta_gjr <-
    rnorm(length(beta_g[[j]][,r]), mean_beta, sqrt(diag(var)))
  # if (j > 1 && r > 1) {
  #   beta_gjr[1, r] <-
  #     rtruncnorm(
  #       1,
  #       b = beta_g[[j]][1, (r - 1)],
  #       mean = (mean_beta)[1],
  #       sd = sqrt(var)[1]
  #     )
  # }
  return(beta_gjr)
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
  sumabsb <- sapply(beta_g, function(beta_gj) {
    apply(beta_gj, 2, function(beta_gjr) {
      sum(abs(beta_gjr))
    })
  })
  if(rank == 1){
    lambda.g <- sapply(sumabsb,function(sumabsb_j){
      rgamma(rank,
             a.lambda + sapply(beta_g,nrow),
             b.lambda + (phi.g * tau.g) ^ (-.5) * sumabsb_j)
    })
    lambda.g <- t(lambda.g)
  }else{
    lambda.g <-
      apply(sumabsb, 2, function(sumabsb_j) {
        rgamma(rank,
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
#' @return an update for omega
#' @keywords internal
BTRR_GGM_draw_omega <- function(beta_g, tau.g, phi.g, lambda.g) {
  omega_g <- mapply(function(beta_gj, lambda_gj) {
    omega_gj <- mapply(function(beta_gjr, phi_gr, lambda_gjr) {
      chi <- beta_gjr ^ 2 / (tau.g * phi_gr)
      omega_gjr <- GIGrvg::rgig(nrow(beta_gj), 0.5, chi, lambda_gjr)
    }, beta_gjr = split(beta_gj, col(beta_gj)), phi_gr = phi.g,
    lambda_gjr = lambda_gj, SIMPLIFY = "array")
  }, beta_gj = beta_g, lambda_gj = split(lambda.g, row(lambda.g)),
  SIMPLIFY = FALSE)
  return(omega_g)
}

#' Draw Phi_g under the stick-breaking prior
#'
#' @param beta_g region beta
#' @param W_g region omega
#' @param Xi_g last draw region Xi
#' @param accept_g acceptance count for region g
#' @param cov_Metro_g metropolis covariance for region g
#' @param tau_g region tau
#' @param alpha region alpha
#'
#' @importFrom abind abind
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
    old_Xi_g <- Xi_g
    if (rank == 1) {
      phi.g <- 1
      accept <- 1
      Xi_g <- 1
    } else{
      accept <- accept_g
      new_Xi_g <-
        c(old_Xi_g + cov_Metro_g %*% rnorm(rank - 1))
      while (length(new_Xi_g[new_Xi_g <= 0]) > 0) {
        new_Xi_g <- c(old_Xi_g + cov_Metro_g %*% rnorm(rank - 1))
      }
      new_post_dens <- sum(sapply(seq(rank - 1), function(cr) {
        stick_break_log_posterior(
          Xi_g = new_Xi_g,
          current_rank =  cr,
          BtWB_g =  BtWB_g,
          p_g = p_g,
          tau_g =  tau_g,
          alpha_g = alpha_g
        )
      }))
      old_post_dens <- sum(sapply(seq(rank - 1), function(cr) {
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

#' Update tau
#'
#' @param a.tau hyperparameter
#' @param b.tau hyperparameter
#' @param p_g region dimension
#' @param phi.g region rank variation weights
#' @param BtWB_g quadratic
#'
#' @return update for tau
#' @keywords internal
BTRR_GGM_draw_tau <- function(a.tau, b.tau, p_g, phi.g, BtWB_g) {
  if(rank == 1){
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
#' @param rank number of ranks in the model
#' @param n_iter number of iterations for the MCMC
#'
#' @return preallocated result list object
#' @keywords internal
BTRR_GGM_empty_results <- function(p, rank, n_iter) {
  G <- length(p)
  check_D <- sapply(p,length)
  if(any(check_D !=  check_D[1]))
    stop("All of the regions should have the same dimension length.")
  D <- unique(check_D)
  out <- list(
    betas = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(D), function(j) {
          matrix(NA, p[[g]][j], rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    W = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(D), function(j) {
          matrix(NA, p[[g]][j], rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    lambda = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        sapply(seq(D), function(j) {
          matrix(NA, D, rank)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F),
    Phi = sapply(seq(n_iter), function(s) {
      sapply(seq(G), function(g) {
        rep(NA,rank)
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
#' @return A list of the elements from the most recent MCMC draw
#' @keywords internal
BTRR_GGM_extract_last_iter <- function(result_file) {
  result_list <- readRDS(result_file)
  # Backward compatibility fix
  if("B" %in% names(result_list)) result_list$betas <- result_list$B
  out <- list()
  out$S <- length(results_list$betas)
  out$d <- results_list$d[,,S]
  out$betas <- results_list$betas[[S]]
  G <- length(out$betas)
  if(G > 1){
    out$Sig <- results_list$Sig[,,S]
    out$zeta <- results_list$zeta[S]
    out$Upsilon <- matrix(1, G, G) ## Suggested by Wang
    diag(out$Upsilon) <- 0 ## Suggested by Wang
    for(g in seq(G)){
      if(g < G){
        for(gg in (g+1):G){
          out$Upsilon[gg,g] <- out$Upsilon[g,gg] <- 1/rinvgauss(1,sqrt(out$zeta^2 / out$Sig[g,gg]^2),out$zeta^2)
        }
      }
    }
  }
  out$tau <- results_list$tau[,S]
  # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
  out$alpha.g <-
    seq(rank ^ (-D), rank ^ (-.1), length.out = 10)
  out$lambda <- results_list$lambda[[S]]
  out$W <- results_list$W[[S]]
  out$Phi <- results_list$Phi[[S]]
  out$sig2y <- results_list$sig2y[S]
  if(rank > 1){
    out$Xi <- sapply(out$Phi,function(each_rank){
      out <- numeric(rank - 1)
      out[1] <- each_rank[1]
      for(r in seq(2,rank - 1)){
        out[r] <- each_rank[r] / (prod(head((1 - each_rank),r-1)))
      }
      return(out)
    })
  }else{
    out$Xi <- matrix(1,1,G)
  }
  out$cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(rank - 1),simplify = FALSE)
  return(out)
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
    require(abind)
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

# OLD FUNCTIONS ----
make_region_array <- function(y){
  lapply(sapply(1:G, function(each_region) {
    lapply(y, function(each_subject) {
      each_subject[[each_region]]
    })
  }, simplify = F), function(each_subject_and_region) {
    array(unlist(each_subject_and_region), dim = c(
      dim(each_subject_and_region[[1]]),
      length(each_subject_and_region)
    ))
  })
}

single_region_array <- function(y, region) {
  sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array")
}

std_single_region_array <- function(y,region,Dim,TT,p,n){
  a <- sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array")
  b <- apply(a,Dim+2,function(each_subject){
    y_mean <- apply(each_subject,seq(Dim),mean)
    y_sd <- apply(each_subject,seq(Dim),sd)
    centered_y <- apply(each_subject,Dim+1,function(each_time){
      (each_time - y_mean) / y_sd
    })
    centered_y[is.nan(centered_y)] <- 0
    array(centered_y,dim = dim(each_subject)[seq(Dim+1)])
  })
  cc <- array(b,dim = c(p[[region]],TT,n))
  return(cc)
}

vec_std_single_region_array <- function(y,region,p){
  sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array") %>%
    apply(D+2,function(each_subject){
      y_mean <- apply(each_subject,seq(D),mean)
      y_sd <- apply(each_subject,seq(D),sd)
      centered_y <- apply(each_subject,D+1,function(each_time){
        (each_time - y_mean) / y_sd
      })
      centered_y[is.nan(centered_y)] <- 0
      array(centered_y,dim = c(p[[region]],TT))
    }) %>%
    array(dim = c(prod(p[[region]]),TT,n))
}

#' Extract tensor coefficient elements
#'
#' @param results_B The \code{B} list element in an object output from \code{BTR_Y_x_with_GGM}
#'
#' @return A list of arrays, each with a given region's posterior samples of the tensor coefficient in an array.
#' @export
#'
#' @examples
#' \dontrun{
#' random_seed <- 610
#' set.seed(random_seed)
#' input <-
#'   Y_x_with_GGM_fmri_data(
#'     subjects = 20,
#'     regions = 10,
#'     max_time = 100,
#'     avg_margin_size = c(10, 10, 10),
#'     SNR = 5,
#'     CNR = 1,
#'     conn_regions = 2,
#'     conn_level = 0.9,
#'     block_period = 30,
#'     scenario = NULL
#'   )
#' test_output <- BTR_Y_x_with_GGM(
#' input = input,
#' n.iter = 110,
#' n.burn = 0,
#' ranks = 1,
#' hyperparameters =
#'   data.frame(
#'     a.tau = D - 1,
#'     b.tau = rr^((1 / D) - 1), # These values are from Guhaniyogi et al [2017]
#'     a.lambda = 3,
#'     b.lambda = 3^(1/(2*D)), # These values are from Guhaniyogi et al [2017]
#'     a.zeta = 1,
#'     b.zeta = 0.01,
#'     a.sig = 1,
#'     b.sig = -log(0.95) # These values are from Guhaniyogi et al [2017]
#'   )
#' ,
#' save_after = 1e5,
#' save_llik = TRUE
#' )
#' all_B <- extract_tensor_coefficient(test_output$B)
#' }
extract_tensor_coefficient <- function(results_B) {
  out <- sapply(seq(length(results_B[[1]])),function(each_region) {
    sapply(seq(length(results_B)),function(each_iter) {
      composeParafac(results_B[[each_iter]][[each_region]])
    },simplify = "array")
  },simplify = FALSE)
  return(out)
}

asplit <- function(X,margin) {
  if(length(dim(X)) < 2) stop("This will only work matrices or arrays of dimension 2 or higher.")
  if(margin > length(dim(X))) stop("Dimension of X is less than margin.")
  X_dim <- dim(X)
  split_list <- split(X,slice.index(X,margin))
  output <- lapply(split_list,function(j){
    out <- array(j,dim = X_dim[-margin])
    return(out)
  })
  return(output)
}

lb_B <- function(b, alpha = 0.05) {
  out <- apply(b, seq(length(dim(b)) - 1), quantile, probs = alpha/2)
  return(out)
}

ub_B <- function(b, alpha = 0.05) {
  out <- apply(b, seq(length(dim(b)) - 1), quantile, probs = 1 - alpha/2)
  return(out)
}

CR_contain_zero <- function(b, alpha = 0.05) {
  lb <- lb_B(b,alpha)
  ub <- ub_B(b,alpha)
  out <- as.numeric(lb > 0 | ub < 0)
  return(out)
}

point_estimate_credible <- function(b, alpha = 0.05) {
  significant <- CR_contain_zero(b,alpha)
  med_B <- apply(b,seq(length(dim(b)) - 1),median)
  out <- array(significant,dim = dim(med_B)) * med_B
  return(out)
}

s2m <- function(x,b){
  two_means <- try(kmeans(abs(x),2))
  zero_idx <- which(two_means$cluster == which.min(two_means$centers))
  A <- x[zero_idx]
  if(length(A) < 2) {
    return(1)
    break
  }
  two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
  iterations <- 1
  while(abs(two_centers$centers[1, 1] - two_centers$centers[2, 1]) > b) {
    zero_idx <- which(two_centers$cluster == which.min(two_centers$centers))
    A <- A[zero_idx]
    if(length(A) < 2) break
    two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
    iterations <- iterations + 1
  }
  num_nonzero <- length(x) - length(A)
  return(num_nonzero)
}

s2m_B <- function(B,sigma = NULL){
  if(is.null(sigma)) {
    sigma <- median(apply(B,seq(length(dim(B)) - 1),sd))
  }
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}

s2m_allB <- function(all_B, sigma=NULL) {
  require(abind)
  bound_regions <- abind(all_B, along = 0)
  bound_final_estimates <- s2m_B(bound_regions)
  out <- asplit(bound_final_estimates, margin = 1)
  return(out)
}

#' Make a mask object for use in stitch brain functions
#'
#' @param file_location The relative file path to a binarized mask file with file extension \code{.nii.gz}.
#' @param slice_id The slice along the z-axis that the mask is being applied to.
#'
#' @return A list with three elements \code{Name}, \code{Location}, and \code{Binary}.
#'
make_mask <- function(file_location,slice_id = NULL) {
  if(is.null(slice_id)) stop("This function currently only supports 2D masks. Please specify a slice ID number.")
  require(oro.nifti)
  roi_text <- gsub("_"," ",file_location)
  roi_match <- regexec("(?<first>[[:upper:]][[:lower:]]+) (?<last>[[:upper:]][[:lower:]]+)*",text=roi_text,perl = TRUE)
  roi_name <- trimws(regmatches(roi_text,roi_match)[[1]][1])
  overall_map <- oro.nifti::readNIfTI(file_location)@.Data[,,slice_id+1]
  cut_map <- matrix(0,2,2)
  for(dimm in 1:2){
    bins <- unlist(apply(overall_map,dimm,function(x) which(x == 1)))
    cut_map[dimm,] <- c(min(bins),max(bins))
  }
  cut_roi <- overall_map[seq(cut_map[2,1],cut_map[2,2]),seq(cut_map[1,1],cut_map[1,2])]
  # cut_roi[cut_roi == 0] <- NA
  out <- list(Name = roi_name,
              Location = cut_map,
              Binary = cut_roi)
  return(out)
}

# Sensitivity and specificity funtion
Bayes_sens_spec <- function(B_estimate, true_B) {
  total_voxels <- sum(sapply(true_B,length))
  total_true_pos <- sum(sapply(true_B, function(Bg) {
    sum(Bg != 0)},simplify = T))
  total_true_neg <- sum(sapply(true_B, function(Bg) {
    sum(Bg == 0)}, simplify = T))
  est_true_pos <- sum(mapply(function(est,tru) {
    where_true_pos <- which(tru != 0)
    return(sum(est[where_true_pos] != 0))
  },est = B_estimate, tru = true_B,SIMPLIFY=T))
  est_true_neg <- sum(mapply(function(est,tru) {
    where_true_neg <- which(tru == 0)
    return(sum(est[where_true_neg] == 0))
  },est = B_estimate, tru = true_B,SIMPLIFY=T))
  sensitivity = est_true_pos / total_true_pos
  specificity = est_true_neg / total_true_neg
  return(c(Sensitivity = sensitivity, Specificity = specificity))
}

# DIC function
bayes_tensor_dic <- function(mdl, input) {
  mdl$llik <- mdl$llik[mdl$llik !=0]
  n_burned <- length(mdl$llik) - length(mdl$B)
  D <- ifelse(n_burned == 0,-2*mdl$llik,-2*mdl$llik[-seq(n_burned)])
  if(class(mdl$B[[1]][[1]]) == "list") {
    B_mcmc <- extract_tensor_coefficient(mdl$B)
    posterior_mean_B <- sapply(B_mcmc, function(b_mcmc) {
      out <- apply(b_mcmc,seq(length(dim(b_mcmc)) - 1),mean)
      return(array(out,dim = head(dim(b_mcmc),-1)))
    }, simplify = F)
  } else {
    B_mcmc <- sapply(seq(length(mdl$B[[1]])), function(g) {
      sapply(seq(length(mdl$B)),function(s) {
        mdl$B[[s]][[g]]}, simplify = "array")},simplify = F)
    posterior_mean_B <- sapply(B_mcmc, FUN = apply, 1, mean,simplify=F)
    posterior_mean_B <- mapply(function(b,y) array(b,dim=head(dim(y),-2)),b = posterior_mean_B, y = input$Y, SIMPLIFY = F)
  }
  posterior_mean_d <- apply(mdl$d,1:2,mean,na.rm=T)
  posterior_mean_sig2y <- mean(mdl$sig2y,na.rm=T)
  llik_theta_hat <- sum(mapply(function(y,b,d) {
    BX <- b %o% input$x
    out <- sum(mapply(function(ygi,BXi,dgi) {
      sum(dnorm(c(ygi),mean = c(BXi + dgi),sd = sqrt(posterior_mean_sig2y),log=T))
    },ygi = asplit(y,length(dim(y))), BXi = asplit(BX,length(dim(BX))), dgi = d))
    return(out)
  },y = input$Y, b = posterior_mean_B, d = asplit(posterior_mean_d,1)))
  p_d <- 2*(llik_theta_hat + mean(D)/2)
  dic <- p_d - 2*llik_theta_hat
  return(dic)
}

# Function to get a point estimate of Sigma
s2m_Sigma <- function(Sig) {
  keep_upper_tri <- upper.tri(Sig[,,1])
  off_diag <- sapply(seq(tail(dim(Sig),1)), function(s) {
    Sig[,,s][keep_upper_tri]
  }, simplify = "array")
  final_off_diag <- s2m_B(off_diag)
  final_on_diag <- diag(apply(Sig,1:2,median))
  final_Sig <- diag(0,length(final_on_diag))
  final_Sig[upper.tri(final_Sig)] <- final_off_diag
  final_Sig <- final_Sig + t(final_Sig)
  final_Sig <- final_Sig + diag(final_on_diag)
  return(final_Sig)
}

# Function to get a point estimate of the partial correlation from Sigma
s2m_part <- function(Sig) {
  require(DensParcorr)
  keep_upper_tri <- upper.tri(Sig[,,1])
  part_corr <- sapply(seq(dim(Sig)[3]), function(s) {
    return(prec2part(Sig[,,s]))
  }, simplify="array")
  off_diag <- sapply(seq(tail(dim(part_corr),1)), function(s) {
    part_corr[,,s][keep_upper_tri]
  }, simplify = "array")
  final_off_diag <- s2m_B(off_diag)
  final_on_diag <- diag(apply(part_corr,1:2,median))
  final_part <- diag(0,length(final_on_diag))
  final_part[upper.tri(final_part)] <- final_off_diag
  final_part <- final_part + t(final_part)
  final_part <- final_part + diag(final_on_diag)
  return(final_part)
}

# This is a function for finding the locations of false positives
fp_loc <- function(activation_estimate,activation_truth) {
  # require(tidyverse)
  require(reshape2)
  est <- melt(activation_estimate,value.name = "estimate")
  tru <- melt(activation_truth, value.name = "truth")
  et_df <- suppressMessages(left_join(est,tru))
  fp_df <- filter(et_df, estimate != 0 & truth == 0)
  out <- select(fp_df, -truth, -estimate)
  return(out)
}

# This function finds the true positive locations
tp_loc <- function(activation_truth) {
  # require(tidyverse)
  require(reshape2)
  out <- melt(activation_truth) %>% filter(value != 0) %>% select(-value)
  return(out)
}

# This is a function for determining the distance from a true activation
activation_distances <- function(fp_locations, tp_locations) {
  fp_split <- split(fp_locations, row(fp_locations))
  min_dists <- sapply(fp_split, function(fp) {
    locs <- rbind(fp, tp_locations)
    out <- min(as.matrix(dist(locs))[1,-1])
    return(out)
  })
  return(min_dists)
}

# This is a function that finds the average distance to activation for a model
model_avg_fp_dist <- function(EST_tensors, TRU_tensors) {
  all_dists <- mapply(function(ee,tt) {
    FPL <- fp_loc(ee, tt)
    TPL <- tp_loc(tt)
    ADs <- activation_distances(FPL,TPL)
    return(ADs)
  }, ee = EST_tensors, tt = TRU_tensors, SIMPLIFY = T)
  return(mean(unlist(all_dists)))
}

# This is a function that finds the median distance to activation for a model
model_median_fp_dist <- function(EST_tensors, TRU_tensors) {
  all_dists <- mapply(function(ee,tt) {
    FPL <- fp_loc(ee, tt)
    TPL <- tp_loc(tt)
    ADs <- activation_distances(FPL,TPL)
    return(ADs)
  }, ee = EST_tensors, tt = TRU_tensors, SIMPLIFY = T)
  return(median(unlist(all_dists)))
}

