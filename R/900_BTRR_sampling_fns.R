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
BTRR_draw_Xi <- function(Xi, cov_Metro, betas, W, tau, a) {
  proposal <- Xi + rnorm(length(Xi))%*%chol(cov_Metro)
  while(any(proposal < 0) | any(proposal > 1)) {
    proposal <- Xi + rnorm(length(Xi))%*%chol(cov_Metro)
  }
  rank <- ncol(betas[[1]])
  ld_Xi <- sum(sapply(Xi,function(xi) sum(dbeta(xi,1,a,log = TRUE))))
  ld_proposal <- sum(sapply(proposal,function(xi) sum(dbeta(xi,1,a,log = TRUE))))
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
BTRR_draw_tau <- function(a.tau,b.tau,betas,W,Phi) {
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
#' @importFrom truncnorm rtruncnorm
#'
#' @return A scalar
#' @keywords internal
BTRR_draw_k <- function(Y,x,betas,sigma_epsilon_sq) {
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
BTRR_draw_sigma_epsilon_sq <- function(a_epsilon,b_epsilon,Y,x,betas,k) {
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
BTRR_draw_beta <- function(Y,x,betas,Sigma_AR,tau,Phi,W,j,r) {
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
composeParafac2 <- function(bb){
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
composeParafac <- function(bb) {
  Rank <- ncol(bb[[1]])
  core_tensor <- diag(1,Rank)
  out <- Reduce(`%x%`, rev(bb)) %*% c(core_tensor)
  return(array(out, dim = sapply(bb,nrow)))
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
