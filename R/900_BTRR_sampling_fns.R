#' Compose a tensor from its CANDECOMP/PARAFAC (CP) decomposition
#'
#' This function takes a list of length D containing all of
#' the components of the CP decomposition and returning a D-dimensional
#' tensor.
#'
#' @param bb A list of length D in which each element is a p_d by R matrix
#'
#' @return A single array-class tensor. In two-dimensions, this will be returned
#'   as a matrix.
#' @export
#'
composeParafac <- function(bb){
  DD <- length(bb)
  pp <- lapply(bb,nrow)
  RR <- lapply(bb,ncol)
  stopifnot(all(unlist(RR) == unlist(RR)[1]))
  RRR <- unique(unlist(RR))
  cp_summands <- sapply(1:RRR,function(r){
    rank_betas <- lapply(bb,function(b){b[,r]})
    Reduce("%o%",rank_betas)
  },simplify = "array")
  apply(cp_summands,1:DD,sum)
}

#' Matricization of a tensor
#'
#' @param tens An array of dimension \eqn{D}
#' @param k An integer less than or equal to \eqn{D}, which will become the
#'   first index in the matricized array
#'
#' @return A matrix
#' @keywords internal
mode_k_matriz <- function(tens,k){
  if(k > length(dim(tens))) stop("You cannot have a mode greater than the dimension of the tensor!")
  t(apply(tens,k,c))
}

#' Calculate stick values
#'
#' @param breaks A vector of length \eqn{n-1} of break proportions
#'
#' @return A vector of length \eqn{n} of break locations between 0 and 1
#' @keywords internal
stick_values <- function(breaks){
  if(!is.numeric(breaks)) stop("You need to input a numeric variable.")
  out <- numeric(length(breaks) + 1)
  for(i in seq(length(breaks))){
    out[i] <- breaks[i]*prod(1 - breaks[seq(length(breaks)) < i])
  }
  out[length(breaks) + 1] <- 1 - sum(out[seq(length(breaks))])
  return(out)
}

#' The log of the posterior probability for the stick-breaking prior structure
#'
#' @param Xi A vector of length \eqn{R - 1}
#' @param betas,W Lists of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param tau A scalar value for the global variance
#' @param a A scalar value for the value of variance weight
#'
#' @return An integer value
#' @keywords internal
stick_break_log_posterior <- function(Xi,betas,W,tau,a) {
  R <- length(Xi) + 1
  part_a <- -sum((sum(sapply(betas,nrow)) / 2) * log(Xi))
  part_b <- sapply(seq(length(Xi)),function(r) {(a -(R - r)*
      sum(sapply(betas,nrow)) / 2 - 1) * log(1-Xi[r])})
  Cr <- mapply(function(b_j,W_j){
    mapply(function(b_jr,W_jr){
      crossprod(b_jr,diag(1/W_jr)) %*% b_jr
    },b_jr = split(b_j,col(b_j)), W_jr = split(W_j,col(W_j)))
  },b_j = betas, W_j = W)
  part_c <- -sum((1/Xi)*head(rowSums(Cr),-1))/tau
  part_d <- 0
  if(R > 2){
    part_d <- -sum(sapply(seq(2,R-1),
                         function(k){
                           (1/Xi[k]*prod(1 - Xi[seq(k-1)]))*
                             rowSums(Cr)[k]/tau}))
  }
  part_e <- -(1/(prod(1 - Xi)*tau))*rowSums(Cr)[R]
  out <- sum(c(part_a,part_b,part_c,part_d,part_e))
  return(out)
}

#' Draw Xi
#'
#' Perform a Metropolis-Hastings draw of the Xi parameters.
#'
#' @param Xi The previous iteration values for Xi
#' @param cov_Metro The Metropolis proposal covariance
#' @param betas,W Lists of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param tau A scalar value for the global variance
#' @param a A scalar value for the value of variance weight
#'
#' @return A vector of length \eqn{R - 1}
#' @keywords internal
draw_Xi <- function(Xi, cov_Metro, betas, W, tau, a) {
  proposal <- Xi + rnorm(length(Xi))%*%chol(cov_Metro)
  while(any(proposal < 0) | any(proposal > 1)) {
    proposal <- Xi + rnorm(length(Xi))%*%chol(cov_Metro)
  }
  ld_Xi <- stick_break_log_posterior(Xi,betas,W,tau,a)
  ld_proposal <- stick_break_log_posterior(proposal,betas,W,tau,a)
  out <- Xi
  if(runif(1) < exp(ld_proposal - ld_Xi)){
    out <- proposal
  }
  return(c(out))
}

#' Draw tau
#'
#' Draw a value of tau using its posterior full conditional distribution.
#'
#' @param a.tau,b.tau Scalar values for the hyperparameters
#' @param betas,W Lists of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param Phi A vector whose elements sum to 1
#'
#' @return A scalar
#' @keywords internal
draw_tau <- function(a.tau,b.tau,betas,W,Phi) {
  D <- length(betas)
  R <- ncol(betas[[1]])
  a_post <- a.tau + D*R/2
  quad_term <- mapply(function(b_j,W_j){
    mapply(function(b_jr,W_jr,p_r){
      crossprod(b_jr,diag(1/W_jr)) %*% b_jr
    },b_jr = split(b_j,col(b_j)), W_jr = split(W_j,col(W_j)), p = Phi)
  },b_j = betas, W_j = W)
  b_post <- b.tau + sum(quad_term)/2
  tau_out <- 1/rgamma(1,a_post,b_post)
  return(tau_out)
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

#' Draw k
#'
#' Draw k from its posterior full conditional distribution.
#'
#' @param Y Array of the response values
#' @param x Matrix of the covariates
#' @param betas List of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param sigma_epsilon_sq Scalar observational variance
#'
#' @return A scalar
#' @keywords internal
draw_k <- function(Y,x,betas,sigma_epsilon_sq) {
  requireNamespace("truncnorm")
  D <- length(betas)
  time_T <- dim(Y)[D + 1]
  B <- composeParafac(betas)
  epsilon <- apply(Y - (B %o% x),D + 1,identity)
  sum_epsilon_t_minus_1_sq <- sum(epsilon[,seq(time_T - 1)]^2)
  sum_epsilon_t_epsilon_t_minus_1 <- sum(epsilon[,seq(2,time_T)]*epsilon[,seq(time_T - 1)])
  out <- truncnorm::rtruncnorm(1,a = -1, b = 1,
                    mean = sum_epsilon_t_epsilon_t_minus_1 /
                      sum_epsilon_t_minus_1_sq,
                    sd = sqrt(sigma_epsilon_sq / sum_epsilon_t_minus_1_sq))
  return(out)
}

#' Draw sigma_epsilon_sq
#'
#' Draw sigma_epsilon_sq from its posterior full conditional distribution.
#'
#' @param a_epsilon,b_epsilon Scalar values of the hyperparameters
#' @param Y Array of the response values
#' @param x Matrix of the covariates
#' @param betas List of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param k Scalar for the AR(1) parameter
#'
#' @return A scalar
#' @keywords internal
draw_sigma_epsilon_sq <- function(a_epsilon,b_epsilon,Y,x,betas,k) {
  D <- length(betas)
  time_T <- dim(Y)[D + 1]
  B <- composeParafac(betas)
  epsilon <- apply(Y - (B %o% x),D + 1,identity)
  epsilon_sum_sq_diff <- sum((epsilon[,seq(2,time_T)] - k*epsilon[,seq(time_T - 1)])^2)
  out <- 1/rgamma(1,a_epsilon + prod(dim(Y))/2, b_epsilon + epsilon_sum_sq_diff/2)
  return(out)
}

#' Draw betas
#'
#' Draw beta variables from their posterior full conditional distributions. As a
#'   note, this is done with a loop in the code for the Markov Chain Monte
#'   Carlo.
#'
#' @param Y Array of the response values
#' @param x Matrix of the covariates
#' @param betas,W Lists of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
#'   matrices
#' @param Sigma_AR A covariance matrix based off of the observational variance \
#'   and the covariance structure.
#' @param tau A scalar for the global variance
#' @param Phi A vector for the rank-weighting of the variance
#' @param j Dimension index
#' @param r Rank index
#'
#' @return A list of length \eqn{D} made up of matrices with dimensions
#'   \eqn{p_j} by \eqn{R}.
#' @keywords internal
draw_beta <- function(Y,x,betas,Sigma_AR,tau,Phi,W,j,r) {
  xtx <- crossprod(x,chol2inv(chol(Sigma_AR))) %*% x
  B_r <- sapply(betas,function(b,r) b[,r,drop = FALSE],r=r,simplify = FALSE)
  if(length(betas) < 3){
    var_beta <- 1 / (sum(B_r[[-j]]^2) * c(xtx) + (1/(tau*Phi[r]) * (1/W[[j]][,r])))
  }else{
    var_beta <- 1 / (sum(Reduce(`%o%`,B_r[-j])^2) * c(xtx) +
      (1/(tau*Phi[r]) * (1/W[[j]][,r])))
  }
  if(ncol(betas[[1]]) == 1){
    Y_til <- Y
  }else{
    B_not_r <- composeParafac(sapply(betas,function(b,r) b[,-r,drop = FALSE], r=r,
                                     simplify = FALSE))
    Y_til <- Y - B_not_r %o% x
  }
  beta_almost_mean <- apply(Y_til,j, function(Y_til_j){
    Y_Sig_x <- apply(Y_til_j,seq(length(betas) - 1),function(each_Y_t){
      crossprod(each_Y_t,chol2inv(chol(Sigma_AR)))%*%x
    })
    sum(composeParafac(B_r[-j]) * Y_Sig_x)
  })
  beta_mean <- var_beta * beta_almost_mean
  out <- rnorm(length(betas[[j]][,r]),beta_mean,sqrt(var_beta))
  return(c(out))
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
post_hoc_AR1_llik <- function(input_data,results_obj) {
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
