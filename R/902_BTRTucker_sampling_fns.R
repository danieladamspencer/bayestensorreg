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
#' @return updates the betas object at index j
#' @keywords internal
BTRT_draw_Bj <- function(y_til,X,betas,tau,W,sig_y2,G,j) {
  prior_covar <- tau * c(W[[j]])
  B_notj <- Reduce(`%x%`,rev(betas[-j])) %*% apply(G,j,identity)
  Delta <- apply(X,length(dim(X)),function(X_i) {
    t(apply(X_i,j,identity)) %*% B_notj
  })
  covar_Bj <- chol2inv(chol(diag(1/prior_covar) + tcrossprod(Delta) / sig_y2,tol = 1e-18))
  mean_Bj <- covar_Bj %*% (Delta %*% y_til / sig_y2)
  out <- c(mean_Bj) + rnorm(length(mean_Bj)) %*% chol(covar_Bj)
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
#' @return Update to the core tensor
#' @keywords internal
BTRT_draw_G <- function(y_til, betas, X, z, V, sig_y2) {
  prior_covar <- z * c(V)
  Delta <- apply(X,length(dim(X)),function(X_i) {
    crossprod(Reduce(`%x%`,rev(betas)), c(X_i))
  })
  covar_G <- chol2inv(chol(diag(1/prior_covar) + tcrossprod(Delta) / sig_y2, tol = 1e-20))
  mean_G <- covar_G %*% (Delta %*% y_til / sig_y2)
  out <- c(mean_G) + rnorm(length(mean_G)) %*% chol(covar_G)
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
  out <- GIGrvg::rgig(1,a.z - length(G)/2,sum(G^2 / V), 2*b.z)
  return(out)
}
