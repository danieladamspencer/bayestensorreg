#' Construct tensor coefficient from results
#'
#' As the \code{BTRR_single_subject} function saves posterior samples of the
#'   tensor decomposition components, the tensor coefficient is composed from
#'   results of class \code{BTRR_result}.
#'
#' @param btrr_result A list object with class \code{BTRR_result}
#'
#' @return An array in which the last index corresponds to the \eqn{i}th draw
#'   from the posterior distribution.
#' @export
BTRR_extract_B <- function(btrr_result) {
  out <- sapply(btrr_result$betas, composeParafac, simplify = "array")
  return(out)
}

#' Calculate the log-likelihood
#'
#' Calculate the log-likelihood for the Bayesian tensor response regression
#'   model. This can be useful when performing diagnostics. This is excluded
#'   from the Markov Chain function for performance gains.
#'
#' @param input_data A list of class \code{TRR_data} containing (at least) a
#'   response array \code{Y} and a design matrix \code{x}.
#' @param results_obj A list of class \code{BTRR_result} containing the draws
#'   from the posterior distribution from \code{\link{BTRR_single_subject}}.
#'
#' @return A vector with length equal to the number of samples from the
#'   posterior distribution containing the log-likelihood values for each sample.
#' @export
BTRR_post_hoc_AR1_llik <- function(input_data,results_obj) {
  if(!requireNamespace("mvtnorm")) stop("You will need the mvtnorm package.")
  S <- length(results_obj$B)
  each_llik <- sapply(seq(S), function(s) {
    XB <- composeParafac(results_obj$B[[s]]) %o% input_data$x
    Sigma_AR <- toeplitz(results_obj$k[s]^seq(0,dim(input_data$Y)[length(dim(input_data$Y)) - 1] - 1)) * (results_obj$sigma_epsilon_sq[s] / (1 - results_obj$k[s]^2))
    Y_til <- input_data$Y - XB
    out <- sum(apply(Y_til,seq(length(dim(input_data$Y)) - 2), function(y_til) {
      dmvnorm(c(y_til),sigma = Sigma_AR, log = TRUE)
    }))
  })
  return(each_llik)
}


#' Compose all tensors with Tucker decomposition
#'
#' Uses the tensor decomposition from Tucker (1927) to compose a tensor
#'   coefficient draw from a \code{BTRT_result} object.
#'
#' @param btrt_object Output object from \code{BTRTucker}
#'
#' @return An array of samples from the posterior distribution of the tensor
#'   coefficient. The last dimension corresponds to the draw number
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- BTRTucker(input)
#' all_B <- BTRT_all_B(results)
#' }
BTRT_all_B <- function(btrt_object) {
  all_B <- sapply(seq(length(btrt_object$betas)), function(iter_idx) {
    composeTuckerCore(btrt_object$betas[[iter_idx]],btrt_object$G[[iter_idx]])
  }, simplify = "array")
  return(all_B)
}

#' Produce a final tensor coefficient estimate
#'
#' This uses the sequential 2-means post-hoc variable-selection technique from
#'   Li and Pati (2018) to produce a point estimate of the tensor coefficient
#'   from all of posterior samples of the tensor coefficient.
#'
#' @param all_B An array of tensor coefficient posterior draws in which the
#'   final dimension corresponds to the posterior sample index.
#'
#' @return An array as a final estimate of the tensor coefficient
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- BTRTucker(input)
#' all_B <- BTRT_all_B(results)
#' estimate_B <- BTRT_estimate_B(all_B)
#' }
BTRT_estimate_B <- function(all_B) {
  B_sd <- median(apply(all_B,seq(length(dim(all_B)) - 1),sd))
  final_B <- s2m_B(all_B, sigma = B_sd)
  return(final_B)
}

#' Produce a final estimate of the tensor coefficient
#'
#' This is a helper function that combines \code{BTRT_all_B},
#'   \code{BTRT_estimate_B}, \code{composeTuckerCore}, \code{s2m_B}, and
#'   \code{s2m}.
#'
#' @param btrt_object An object of class \code{BTRT_result} output from
#'   \code{BTRTucker}
#'
#' @return An array estimate of the tensor coefficient
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' result <- BTRTucker(input)
#' final_B <- BTRT_final_B(result)
#' }
BTRT_final_B <- function(btrt_object) {
  all_B <- BTRT_all_B(btrt_object)
  out <- s2m_B(all_B)
  return(out)
}

#' Compose a tensor based on the Tucker Decomposition
#'
#' @param bb a list of tensor decomposition matrices in which the number of
#'   columns corresponds to the dimension margins
#' @param gg an array in which the length of each dimension margin corresponds
#'   to the number of ranks to decompose each dimension
#'
#' @return an array (tensor) composed using the Tucker tensor decomposition
#' @export
composeTuckerCore <- function(bb,gg) {
  out <- Reduce(`%x%`,rev(bb)) %*% c(gg)
  return(array(out,dim = sapply(bb,nrow)))
}

#' Deviance Information Criterion
#'
#' Calculate the Deviance Information Criterion (DIC) through an approximation
#'   based solely on the values of the log-likelihood at each iteration of the
#'   Markov Chain Monte Carlo results.
#'
#' @param log_likelihood A vector of the log-likelihood values
#' @param burn_in A scalar indicating how many of the values should be removed from the beginning of the \code{log_likelihood} vector
#'
#' @return The DIC approximation value (a scalar)
#' @export
DIC <- function(log_likelihood,burn_in = 0) {
  if(burn_in == 0) {
    D <- -2*log_likelihood
  } else {
    D <- -2*log_likelihood[-seq(burn_in)]
  }
  p_d <- var(D)/2
  dic <- p_d + mean(D)
  return(dic)
}

#' Sequential 2-means number of nonzero-values parameters for an iteration
#'
#' @param x A vector of parameters from a single draw from a posterior
#'   distribution
#' @param b The maximum distance at which two clusters of parameters are
#'   considered equal
#'
#' @return A scalar representing the number of nonzero-valued parameters in a
#'   particular draw from a posterior distribution
#' @export
#'
#' @examples
#' post_sample <- numeric(100)
#' for(k in 1:100) {
#'   g <- runif(1)
#'   if(g < 0.8) {
#'     post_sample[k] <- rnorm(1,0,0.2)
#'   } else {
#'     post_sample[k] <- rnorm(1,0,0.2)
#'   }
#' }
#' how_many_nonzero <- s2m(post_sample, b = sd(post_sample))
s2m <- function(x,b){
  two_means <- kmeans(abs(x),2)
  zero_idx <- which(two_means$cluster == which.min(two_means$centers))
  A <- x[zero_idx]
  two_centers <- try(kmeans(abs(A),2,algorithm=c("Lloyd")))
  if(class(two_centers) != "try-error") {
    iterations <- 1
    while(abs(two_centers$centers[1, 1] - two_centers$centers[2, 1]) > b) {
      zero_idx <- which(two_centers$cluster == which.min(two_centers$centers))
      A <- A[zero_idx]
      two_centers <- try(kmeans(abs(A),2,algorithm=c("Lloyd")))
      if(class(two_centers) == "try-error") break
      iterations <- iterations + 1
    }
  }
  num_nonzero <- length(x) - length(A)
  return(num_nonzero)
}

#' Perform sequential 2-means variable selection on posterior draws for B
#'
#' This is a wrapper function to perform the sequential 2-means variable
#'   selection method [@li2017variable] on the reconstructed tensor from the
#'   posterior draws of betas.
#'
#' @param B The reconstructed tensor from the posterior draws of betas. This
#'   should have dimension \eqn{D + 1}.
#' @param sigma The value for the distance at which two cluster centers are
#'   determined to be the same. This value defaults to the median value for
#'   the posterior standard deviation of the voxel-wise draws in \code{B}.
#'
#' @return An array of dimension \eqn{D + 1}.
#' @export
#'
#' @examples
#' # In this case, the value for D is 2
#' true_B <- rbind(
#'   c(0,0,0,0,0),
#'   c(0,0,1,0,0),
#'   c(0,1,1,1,0),
#'   c(0,0,1,0,0),
#'   c(0,0,0,0,0)
#' )
#' sim_B <- array(NA, dim = c(5,5,500)) # Simulating 500 draws
#' for(s in 1:500) {
#'   sim_B[,,s] <- rnorm(25,mean = c(true_B), sd = 0.3)
#' }
#' final_estimate <- s2m_B(sim_B)
s2m_B <- function(B,sigma = NULL){
  if(is.null(sigma)) sigma <- median(apply(B,seq(length(dim(B)) - 1),sd))
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}
