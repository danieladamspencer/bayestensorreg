#' Obtain the log-likelihood for the tensor regression linear model
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param B tensor coefficient
#' @param gam vector coefficients
#'
#' @importFrom stats dnorm
#'
#' @return scalar
#' @keywords internal
ftr_log_likelihood <- function(input,B,gam) {
  if(!is.null(input$eta)) {
    gam_eta <- c(tcrossprod(gam,input$eta))
  } else {gam_eta <- 0}
  XB <- apply(input$X,length(dim(input$X)),function(x){
    crossprod(c(x),c(B))
  })
  y_res <- input$y - gam_eta - XB
  out <- sum(dnorm(y_res, sd = sd(y_res),log = TRUE))
  return(out)
}

#' Obtain the log-likelihood for the tensor regression probit model
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param B tensor coefficient
#' @param gam vector coefficients
#'
#' @importFrom stats dbinom
#'
#' @return scalar
#' @keywords internal
ftr_probit_log_likelihood <- function(input,B,gam) {
  Lam <- c(tcrossprod(gam,input$eta)) -
    apply(input$X,length(dim(input$X)),function(x){
      crossprod(c(x),c(B))
    })
  PI <- exp(Lam)/(1 + exp(Lam))
  out <- sum(dbinom(input$y,size = 1,prob = PI, log = TRUE))
  return(out)
}

#' Logit transformation
#'
#' @param p vector
#'
#' @return vector
#' @keywords internal
logit <- function(p) {
  log(p/(1-p))
}

#' Inverse logit transformation
#'
#' @param x vector
#'
#' @return vector
#' @keywords internal
expit <- function(x) {
  exp(x)/(1 + exp(x))
}

#' Compose a tensor from its PARAFAC decomposition
#'
#' @param betas a list, where the first level indexes dimension, and the columns
#'   of the matrices within each list represent the rank
#'
#' @return an array
#' @keywords internal
compose_parafac <- function(betas) {
  if(length(unique(sapply(betas,ncol))) != 1) stop("betas should be a list of D matrices. All of these matrices should have the same number of columns.")
  betas_by_rank <- sapply(seq(ncol(betas[[1]])),function(r){
    sapply(betas,function(b_j){
      b_j[,r]
    },simplify = FALSE)
  },simplify = FALSE)
  B_summands <- sapply(betas_by_rank,Reduce,
                       f = `%o%`,simplify = "array")
  out <- apply(B_summands,seq(length(dim(B_summands)) - 1),sum)
  return(out)
}

#' Khatri Rao product
#'
#' @param A matrix
#' @param B matrix
#'
#' @return Khatri Rao matrix product
#' @keywords internal
khatri_rao <- function(A,B = NULL) {
  if(is.null(B)) B <- matrix(1,1,ncol(A))
  if(ncol(A) != ncol(B)) stop("A and B should be matrices with the same number of columns.")
  out <- mapply(function(a,b) a %x% b, a = split(A,col(A)), b = split(B,col(B)))
  return(out)
}

#' Compose a tensor from its Tucker decomposition
#'
#' @param bb list
#' @param gg array
#'
#' @return array
#' @keywords internal
compose_tucker_ftr <- function(bb,gg) {
  # rr <- sapply(bb,ncol)
  # if(all(rr == 1)){
  #   output <- Reduce(outer,sapply(bb,c,simplify = FALSE))
  # }else{
  #   rank_combos <- expand.grid(sapply(rr,seq))
  #   if(length(unique(rr)) == 1) rank_combos <- as.matrix(t(rank_combos))
  #   tucker_summands <- mapply(function(which_ranks,g) {
  #     Reduce(outer,mapply(function(b,r){b[,r]},b = bb, r = which_ranks,SIMPLIFY = FALSE)) * g
  #   },which_ranks = split(as.matrix(rank_combos),row(rank_combos)),g = c(gg),SIMPLIFY = FALSE)
  #   output <- Reduce(`+`,tucker_summands)
  # }
  all_dims <- sapply(bb, nrow)
  vec_out <- compose_tucker_ftr_vec(bb,gg)
  output <- array(vec_out, dim = all_dims)
  return(output)
}

#' Compose a tensor from the Tucker decomposition and vectorize
#'
#' @param bb list
#' @param gg array
#'
#' @return vector
#' @keywords internal
compose_tucker_ftr_vec <- function(bb,gg) {
  out <- Reduce(`%x%`,rev(bb)) %*% c(gg)
  return(out)
}

#' Get the BIC for an FTR model fit
#'
#' @param ftr_object The result from a call to FTR_CP or FTRTucker
#' @param input list with elements \code{y} and \code{X}
#'
#' @return scalar
#' @keywords internal
ftr_model_BIC <- function(ftr_object, input) {
  if(!is.null(ftr_object$G)) {G_parms <- prod(dim(ftr_object$G))} else {G_parms <- 0}
  k <- sum(sapply(ftr_object$betas, function(be) prod(dim(be)))) + length(ftr_object$gam) + G_parms
  n <- length(input$y)
  out <- k*log(n) - 2*ftr_object$llik
  return(out)
}
