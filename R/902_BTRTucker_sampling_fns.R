#' Perform mode-k matricization on an array
#'
#' @param x an array
#' @param k an integer or vector of integers
#'
#' @return a matrix with dimensions \code{prod(dim(x)[k]) x prod(dim(X)[-k])}
#' @export
#'
#' @examples
#' set.seed(47408)
#' p <- rep(50,3)
#' N <- 20
#' A <- array(rnorm(prod(p,N)),dim = c(p,N))
#' Ak <- kFold(A,2)
kFold = function(x,k) {
  return(matrix(c(aperm(x,c(k,seq(length(dim(x)))[-k]))),prod(dim(x)[k]),prod(dim(x)[-k])))
}

#' Draw from PFC for betas
#'
#' @param y_til response minus etagam
#' @param X covariate tensor
#' @param betas all betas, for betas from other dimensions
#' @param tau global variance
#' @param W local variance
#' @param sig_y2 observational variance
#' @param G core tensor
#' @param j dimension index
#'
#' @importFrom Matrix chol solve
#'
#' @return updates the betas object at index j
#' @keywords internal
BTRT_draw_Bj <- function(y_til,X,betas,tau,W,sig_y2,G,j) {
  prior_covar <- tau * c(W[[j]])
  N <- length(y_til)
  B_notj <- as.matrix(Reduce(`%x%`,rev(betas[-j])) %*% apply(G,j,identity))
  Delta <- kFold(X,c(j,length(dim(X)))) %*% B_notj |>
    array(c(dim(X)[j],N,ncol(B_notj))) |>
    aperm(c(1,3,2)) |>
    c() |>
    matrix(prod(dim(X)[j],ncol(B_notj)),tail(dim(X),1))
  prec_Bj <- diag(1/prior_covar) + tcrossprod(Delta) / sig_y2
  chol_prec_Bj <- Matrix::chol(prec_Bj)
  mean_Bj <- Matrix::solve(prec_Bj, (Delta %*% y_til / sig_y2))
  Z <- rnorm(length(mean_Bj))
  out <- backsolve(chol_prec_Bj,Z)
  out <- out + mean_Bj
  return(matrix(c(out),nrow(betas[[j]]),ncol(betas[[j]])))
}

#' Draw gamma from PFC
#'
#' @param eta nuisance input
#' @param Sig_0 prior covariance
#' @param mu_gam prior mean
#' @param y_til response minus BX
#' @param sig_y2 observational variance
#'
#' @return gamma update
#' @keywords internal
BTRT_draw_gam <- function(eta,Sig_0,mu_gam,y_til,sig_y2) {
  Sig_gam <- solve(solve(Sig_0) + crossprod(eta)/sig_y2)
  mu_gam <- Sig_gam %*% t(mu_gam%*%solve(Sig_0) + y_til%*%eta / sig_y2)
  out <- mu_gam + chol(Sig_gam) %*% rnorm(ncol(eta))
  return(c(out))
}

#' Draw core tensor from PFC
#'
#' @param y_til response minus gam_eta
#' @param betas tensor decomposition matrices
#' @param X tensor covariate
#' @param z global variance
#' @param V local variance
#' @param sig_y2 observational variance
#'
#' @importFrom Matrix chol solve
#'
#' @return Update to the core tensor
#' @keywords internal
BTRT_draw_G <- function(y_til, betas, X, z, V, sig_y2) {
  prior_covar <- z * c(V)
  B_noG <- Reduce(`%x%`,rev(betas))
  Delta <- crossprod(B_noG,matrix(c(X),ncol = tail(dim(X),1)))
  prec_G <- diag(1/prior_covar) + tcrossprod(Delta) / sig_y2
  chol_prec_G <- Matrix::chol(prec_G)
  mean_G <- Matrix::solve(prec_G,(Delta %*% y_til / sig_y2))
  Z <- rnorm(length(mean_G))
  out <- backsolve(chol_prec_G,Z) + mean_G
  return(array(c(out),dim = sapply(betas,ncol)))
}

#' Update local shrinkage
#'
#' @param a_lam hyperparameter value
#' @param b_lam hyperparameter value
#' @param betas tensor decomposition matrices
#' @param tau global variance
#'
#' @return update to local shrinkage
#' @keywords internal
BTRT_draw_lam <- function(a_lam,b_lam,betas,tau) {
  sapply(betas, function(bj) {
    apply(bj,2,function(bjr) {
      rgamma(1,a_lam + length(bjr),b_lam + sum(abs(bjr)) / sqrt(tau))
    })
  }, simplify = F)
}

#' Draw local variance
#'
#' @param lam local shrinkage parameter
#' @param betas tensor decomposition matrices
#' @param tau global variance
#'
#' @importFrom GIGrvg rgig
#'
#' @return update to local variance
#' @keywords internal
BTRT_draw_omega <- function(lam,betas,tau) {
  mapply(function(bb,ll) {
    R <- length(ll)
    sapply(seq(R),function(r){
      sapply(bb[,r],function(b_l) {
        GIGrvg::rgig(1,1/2, b_l^2 / tau,ll[r]^2)
      })
    })
  },bb = betas, ll = lam,SIMPLIFY = FALSE)
}

#' Update observational variance
#'
#' @param a.sig hyperparameter
#' @param b.sig hyperparameter
#' @param y_til response less effects from covariates
#'
#' @return update to observational variance
#' @keywords internal
BTRT_draw_sig_y2 <- function(a.sig, b.sig, y_til) {
  final_a <- a.sig + length(y_til)/2
  final_b <- b.sig + sum(y_til^2)/2
  out <- 1/rgamma(1,final_a,final_b)
  return(out)
}

#' Update global variance
#'
#' @param a_tau hyperparameter
#' @param b_tau hyperparameter
#' @param betas tensor decomposition matrices
#' @param W local variance
#'
#' @return update to local variance
#' @keywords internal
BTRT_draw_tau <- function(a_tau,b_tau,betas,W) {
  a_post <- a_tau + sum(sapply(betas,function(b) prod(dim(b))))/2
  b_post <- b_tau + sum(unlist(betas)^2 / unlist(W))/2
  out <- 1/rgamma(1,a_post,b_post)
  return(out)
}

#' Update core tensor shrinkage parameter
#'
#' @param a.u hyperparameter
#' @param b.u hyperparameter
#' @param G core tensor
#' @param z global variance
#'
#' @return update to core tensor shrinkage
#' @keywords internal
BTRT_draw_U <- function(a.u,b.u,G,z) {
  U_out <- rgamma(length(G),a.u + 1, b.u + abs(c(G))/sqrt(z))
  out <- array(U_out, dim = dim(G))
  return(out)
}

#' Update to core tensor local variance
#'
#' @param U local shrinkage parameters
#' @param z global shrinkage parameters
#' @param G core tensor
#'
#' @importFrom GIGrvg rgig
#'
#' @return update to core tensor local variance
#' @keywords internal
BTRT_draw_V <- function(U,z,G) {
  V_out <- GIGrvg::rgig(length(U),1/2,c(G)^2 / z, c(U)^2)
  out <- array(V_out, dim = dim(U))
  return(out)
}

#' Update core tensor global variance
#'
#' @param a.z hyperparameter
#' @param b.z hyperparameter
#' @param G core tensor
#' @param V core tensor local variance
#'
#' @importFrom GIGrvg rgig
#'
#' @return update to core tensor global variance
#' @keywords internal
BTRT_draw_z <- function(a.z,b.z,G,V) {
  lambda_gig <- a.z - length(G)/2
  chi_gig <- sum(G^2 / V)
  if(lambda_gig < 0 & chi_gig == 0) chi_gig <- 0.0001
  out <- GIGrvg::rgig(1, lambda_gig, chi_gig, 2*b.z)
  return(out)
}
