#' Frequentist Tensor Regression with the CP decomposition
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param rank The rank of the CP decomposition to be used
#' @param epsilon a value for the stopping rule of the algorithm. Specifically,
#'   this is the upper bound for the differences in the log-likelihood between
#'   two iterations of the algorithm.
#' @param max.iter the maximum number of iterations to use before stopping the
#'   algorithm
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats rnorm
#'
#' @return A list with elements \code{gam} (vector coefficient result),
#'   \code{betas} (tensor decomposition components), \code{B} (the tensor ,
#'   coefficient), and \code{total_time} (time spent to complete the analysis).
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- FTR_CP(input)
#' }
FTR_CP <- function(input, rank = 1, epsilon = 1e-8, max.iter = 1000) {
  start_model <- proc.time()[3]
  if(!is.null(input$eta)) {
    gam_new <- lm(input$y ~ -1 + input$eta)$coefficients
  } else {
    gam_new <- 0
  }
  tensor_P <- head(dim(input$X),-1)
  beta_new <- sapply(tensor_P, function(p_j) matrix(rnorm(p_j*rank,sd = 0.025),p_j,rank),simplify = FALSE)
  llik <- ftr_log_likelihood(input,compose_parafac(beta_new),gam_new)
  new_llik <- llik #- llik*max(epsilon,1000)
  beta_old <- beta_new
  gam_old <- gam_new
  step <- 1
  while((abs(100*(new_llik - llik)/llik) > epsilon & step <= max.iter & (new_llik - llik) > 0) | step == 1) {
    cat("Step",step,"Log-likelihood",new_llik,", % change:",abs((new_llik - llik)/llik)*100,"\n")
    beta_old <- beta_new
    gam_old <- gam_new
    step <- step + 1
    llik <- new_llik
    if(!is.null(input$eta)) {
      step_y <- input$y - c(tcrossprod(gam_new,input$eta))
    } else {step_y <- input$y}
    for(d in seq(length(tensor_P))) {
      # cat(d,"\n")
      step_X <- t(apply(input$X,length(dim(input$X)),function(X_i){
        X_id <- t(apply(X_i,d,identity))
        XB_not_d <- X_id %*% Reduce(khatri_rao,beta_new[-d])
        return(XB_not_d)
      }))
      cv.beta_new <- try(cv.glmnet(step_X,step_y,type.measure = "mse",alpha = 1,
                                   family = "gaussian",intercept = FALSE))
      #cat("Beta class:",class(cv.beta_new),"\n")
      if(class(cv.beta_new) == "cv.glmnet") {
        beta_new[[d]] <- matrix(c(as.matrix(coef(cv.beta_new,s = "lambda.min"))[-1,]),
                                ncol = rank)
      }else{
        cat("LASSO fit for betas FAILED.\n")
        beta_new[[d]] <- matrix(c(lm(step_y ~ -1 + step_X)$coefficients), ncol = rank)
      }
      beta_new[[d]][is.na(beta_new[[d]])] <- 0 # In some cases with MRI data, all scan values are equal to zero.
    }
    B <- compose_parafac(beta_new)
    if(!is.null(input$eta)) {
      gam_new <- lm(input$y - apply(input$X,length(dim(input$X)),function(x){
        crossprod(c(x),c(B))
      }) ~ -1 + input$eta)$coefficients
    }
    new_llik <- ftr_log_likelihood(input,B,gam_new)
  }
  convergence <- step < max.iter
  if(is.null(input$eta)) gam_old <- NULL
  return(
    list(
      gam = gam_old,
      betas = beta_old,
      B = compose_parafac(beta_old),
      llik = llik,
      convergence = convergence,
      total_time = proc.time()[3] - start_model
    )
  )
}
