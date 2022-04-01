#' Draw Phi variable from the posterior full conditional
#'
#' @param n scalar
#' @param alpha scalar
#' @param betas list of matrices
#' @param W list of matrices
#' @param b_tau scalar
#'
#' @return vector or matrix of draws from the posterior full conditional
#'
#' @importFrom GIGrvg rgig
#'
#' @keywords internal
cp_draw_Phi <- function(n,alpha, betas, W, b_tau) {
  if(ncol(betas[[1]]) == 1){
    out <- 1
  }else{
    l <- alpha - sum(sapply(betas,nrow))/2
    ch <- 2*b_tau
    psis <- sapply(seq(ncol(betas[[1]])),function(r) {
      ps <- 2*sum(sapply(seq(length(betas)),function(j) {
        crossprod(betas[[j]][,r],diag(1/W[[j]][,r]))%*%betas[[j]][,r]
      }))
      psi <- GIGrvg::rgig(n,l,ps,ch)
    })
    if(n == 1){
      out <- psis / sum(psis)
    }else{
      out <- apply(psis,2,function(ps) ps / sum(ps))
    }

  }
  return(out)
}

#' Draw tau from the posterior full conditional
#'
#' @param n scalar
#' @param a_tau scalar
#' @param b_tau scalar
#' @param betas list of matrices
#' @param W list of matrices
#' @param Phi vector or matrix
#'
#' @importFrom GIGrvg rgig
#'
#' @return vector
#' @keywords internal
cp_draw_tau <- function(n,a_tau,b_tau,betas,W,Phi) {
  l <- a_tau - sapply(betas,nrow)*ncol(betas[[1]])/2
  ch <- 2*b_tau
  DD <- sapply(seq(ncol(betas[[1]])),function(r) {
    ps <- sum(sapply(seq(length(betas)),function(j) {
      crossprod(betas[[j]][,r],diag(1/W[[j]][,r]))%*%betas[[j]][,r]
    }))
    return(ps)
  })
  out <- GIGrvg::rgig(n,l,2*sum(DD/Phi),ch)
  return(out)
}

#' Draw lam from posterior FC
#'
#' @param a_lam scalar
#' @param b_lam scalar
#' @param betas list of matrices
#' @param Phi vector or matrix
#' @param tau scalar
#'
#' @importFrom stats rgamma
#'
#' @return vector
#' @keywords internal
cp_draw_lam <- function(a_lam,b_lam, betas,Phi,tau) {
  A <- a_lam + sapply(betas,nrow)
  betas_L1 <- t(sapply(betas,function(beta_j){
    apply(beta_j,2,function(beta_jr){
      sum(abs(beta_jr))
    })
  },simplify = "array"))
  if(ncol(betas[[1]]) == 1) betas_L1 <- t(betas_L1)
  out <- sapply(seq(length(betas)),function(j){
    sapply(seq(ncol(betas[[1]])),function(r){
      rgamma(1,A[j],b_lam + betas_L1[j,r]/sqrt(Phi[r]*tau))
    })
  },simplify = FALSE)
  return(out)
}

#' Draw W from posterior FC
#'
#' @param lam vector
#' @param betas list of matrices
#' @param Phi vector or matrix
#' @param tau scalar
#'
#' @importFrom GIGrvg rgig
#'
#' @return list of matrices
#' @keywords internal
cp_draw_W <- function(lam,betas,Phi,tau) {
  out <- mapply(function(lam_j,beta_j){
    as.matrix(mapply(function(beta_jr,phi_r,lam_jr){
      GIGrvg::rgig(length(beta_jr),1/2,lam_jr^2,beta_jr^2/(phi_r*tau))
    },beta_jr = split(beta_j,col(beta_j)),
    phi_r = Phi, lam_jr = lam_j))
  },lam_j = lam, beta_j = betas,SIMPLIFY = FALSE)
  return(out)
}

#' Compose tensor from its parafac decomposition
#'
#' @param betas list of matrices
#'
#' @return array
#' @export
btr_compose_parafac <- function(betas) {
  betas_by_rank <- sapply(seq(ncol(betas[[1]])),function(r) {
    sapply(seq(length(betas)),function(j) {
      betas[[j]][,r]
    },simplify = F)
  },simplify = F)
  out <- sapply(betas_by_rank,function(betas_r){
    Reduce(`%o%`,betas_r)},simplify = FALSE)
  return(Reduce(`+`,out))
}

#' Wrapper to compose tensors from BTR_CP model object
#'
#' @param btr_cp_object results from a call to BTR_CP
#'
#' @return list of arrays
#' @export
btr_cp_all_B <- function(btr_cp_object) {
  out <- sapply(btr_cp_object$betas, btr_compose_parafac, simplify = "array")
  return(out)
}

#' Produce a final estimate of the BTR CP tensor coefficient
#'
#' @param btr_cp_object an object of class \code{BTRCP_result}
#'
#' @return an array of the final tensor estimate found using sequential 2-means
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' result <- BTR_CP(input)
#' final_B <- btr_cp_final_B(result)
#' }
btr_cp_final_B <- function(btr_cp_object) {
  all_B <- btr_cp_all_B(btr_cp_object)
  B_sd <- median(apply(all_B,seq(length(dim(all_B)) - 1),sd))
  final_B <- s2m_B(all_B, sigma = B_sd)
  return(final_B)
}

#' Draw betas from their posterior FCs
#'
#' @param X design array
#' @param betas list of matrices
#' @param sig_y2 scalar
#' @param y vector
#' @param gam vector
#' @param eta design matrix
#' @param W list of matrices
#' @param Phi vector or matrix
#' @param tau scalar
#' @param j (scalar) dimension index
#' @param r (scalar) rank index
#'
#' @importFrom stats rnorm
#' @importFrom Matrix chol solve
#'
#' @return vector
#' @keywords internal
cp_draw_betas <- function(X,betas,sig_y2,y,gam,eta,W,Phi,tau,j,r) {
  # betas_l_r <- sapply(seq(length(betas))[-j],function(l){
  #   betas[[l]][,r]
  # },simplify = FALSE)
  betas_l_r <- sapply(betas[-j], function(x) x[,r], simplify = F)
  B_not_j_r <- c(Reduce(`%o%`,betas_l_r))
  H_jr <- (kFold(X,c(length(dim(X)),j)) %*% c(B_not_j_r)) |>
    c() |> matrix(nrow = tail(dim(X),1))
  # H_jr <- t(apply(X,length(dim(X)),function(X_i){
  #   crossprod(apply(X_i,j,identity), c(B_not_j_r))
  # }))
  if(ncol(betas[[j]]) > 1){
    betas_not_r <- sapply(betas,function(betas_j){betas_j[,-r,drop = F]},simplify = FALSE)
    B_not_r <- btr_compose_parafac(betas_not_r)
    XB_not_r <- c(kFold(X,length(dim(X))) %*% c(B_not_r))
    # XB_not_r <- apply(X,length(dim(X)),function(X_i){crossprod(c(apply(X_i,j,identity)),c(B_not_r))})
  }else{
    XB_not_r <- 0
  }
  y_til <- y - c(eta %*% gam) - XB_not_r
  Sig_inv <- crossprod(H_jr)/sig_y2 + diag(1/W[[j]][,r])/(Phi[r]*tau)
  cholSig_inv <- Matrix::chol(Sig_inv)
  Mu <- Matrix::solve(Sig_inv, (t(H_jr)%*%y_til/sig_y2))
  Z <- rnorm(nrow(betas[[j]]))
  out <- backsolve(cholSig_inv,Z)
  out <- out + c(Mu)
  # Sigma <- chol2inv(chol(crossprod(H_jr)/sig_y2 + diag(1/W[[j]][,r])/(Phi[r]*tau)))
  # Mu <- Sigma%*%t(H_jr)%*%y_til/sig_y2
  # out <- c(Mu) + c(rnorm(nrow(betas[[j]])) %*% chol(Sigma))
  return(out)
}

#' Draw sig_y2 from posterior FC
#'
#' @param nu scalar
#' @param s_02 scalar
#' @param y_til vector
#' @param mu_gam vector
#' @param eta design matrix
#'
#' @importFrom stats rgamma
#'
#' @return scalar
#' @keywords internal
cp_draw_sig_y2 <- function(nu,s_02,y_til,mu_gam,eta) {
  a_sig <- (length(y_til) + nu) / 2
  b_sig <- (nu * s_02 + sum(y_til^2) - crossprod(y_til,eta)%*%mu_gam)
  out <- 1/rgamma(1,a_sig,b_sig)
  return(out)
}

#' Second way to draw sig_y2 from posterior FC
#'
#' @param a.sig scalar
#' @param b.sig scalar
#' @param y.til vector
#'
#' @importFrom stats rgamma
#'
#' @return scalar
#' @keywords internal
cp_draw_sig_y2_two <- function(a.sig, b.sig, y.til) {
  final_a <- a.sig + length(y.til)/2
  final_b <- b.sig + sum(y.til^2)/2
  out <- 1/rgamma(1,final_a,final_b)
  return(out)
}

#' Draw gam from posterior FC
#'
#' @param eta design matrix
#' @param Sig_0 matrix
#' @param y_til vector
#' @param sig_y2 scalar
#'
#' @importFrom stats rnorm
#'
#' @return vector
#' @keywords internal
cp_draw_gam <- function(eta,Sig_0,y_til,sig_y2) {
  Sig_gam <- solve(crossprod(eta) + solve(Sig_0))
  mu <- Sig_gam%*%crossprod(eta,y_til)
  out <- mu + sqrt(sig_y2)*chol(Sig_gam) %*% rnorm(ncol(eta))
  return(out)
}

#' Second way to draw gam from posterior FC
#'
#' @param eta design matrix
#' @param invSig_0 matrix
#' @param mu_gam vector
#' @param y_til vector
#' @param sig_y2 scalar
#'
#' @importFrom stats rnorm
#' @importFrom Matrix chol solve
#'
#' @return vector
#' @keywords internal
cp_draw_gam_two <- function(eta,invSig_0,mu_gam,y_til,sig_y2) {
  invSig_gam <- invSig_0 + crossprod(eta)/sig_y2
  mu_gam <- Matrix::solve(invSig_gam,t(mu_gam%*%invSig_gam + y_til%*%eta / sig_y2))
  cholInvSig_gam <- Matrix::chol(invSig_gam)
  Z <- rnorm(ncol(eta))
  out <- backsolve(cholInvSig_gam,Z)
  out <- out + mu_gam
  # Sig_gam <- solve(solve(Sig_0) + crossprod(eta)/sig_y2)
  # mu_gam <- Sig_gam %*% t(mu_gam%*%solve(Sig_0) + y_til%*%eta / sig_y2)
  # out <- mu_gam + chol(Sig_gam) %*% rnorm(ncol(eta))
  return(c(out))
}

#' Calculate log-likelihood from result object
#'
#' @param btr_object result from BTR_CP
#' @param input input list used to run BTR_CP
#'
#' @importFrom stats dnorm
#'
#' @return vector of log-likelihoods
#' @export
btr_cp_llik <- function(btr_object, input) {
  if(is.null(input$eta)) input$eta <- matrix(0,nrow = length(input$y), ncol = 1)
  all_B <- btr_cp_all_B(btr_object)
  all_B <- asplit(all_B, length(dim(all_B)))
  out <- sapply(seq(length(all_B)), function(s) {
    llik <- sum(dnorm(
      input$y -
        apply(input$X,
              length(dim(input$X)),
              function(X_i)
                crossprod(c(X_i), c(all_B[[s]]))) -
        input$eta %*% as.matrix(btr_object$gam)[s,],
      mean = 0,
      sd = sqrt(btr_object$sig_y2[s]),
      log = TRUE
    ))
    return(llik)
  }, simplify = T)
  return(out)
}
