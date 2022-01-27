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
#' results <- FTR_CP(input$y, input$X, input$eta)
#' }
FTR_CP <- function(input,rank = 1,epsilon = 1e-8) {
  start_model <- proc.time()[3]
  gam_new <- lm(input$y ~ -1 + input$eta)$coefficients
  beta_new <- sapply(head(dim(input$X),-1),function(p_j) matrix(rnorm(p_j*rank,sd = 0.025),p_j,rank),simplify = FALSE)
  llik <- ftr_log_likelihood(input,compose_parafac(beta_new),gam_new)
  new_llik <- llik + max(epsilon,1)
  beta_old <- beta_new
  gam_old <- gam_new
  step <- 1
  while(new_llik - llik > epsilon) {
    cat("Step",step,"Log-likelihood",new_llik,"\n")
    beta_old <- beta_new
    gam_old <- gam_new
    step <- step + 1
    llik <- new_llik
    step_y <- input$y - c(tcrossprod(gam_new,input$eta))
    for(d in seq(length(dim(input$X)) - 1)) {
      cat(d,"\n")
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
        beta_new[[d]] <- matrix(c(lm(step_y ~ -1 + step_X)$coefficients), ncol = rank)
      }
      beta_new[[d]][is.na(beta_new[[d]])] <- 0 # In some cases with MRI data, all scan values are equal to zero.
    }
    B <- compose_parafac(beta_new)
    gam_new <- lm(input$y - apply(input$X,length(dim(input$X)),function(x){
      crossprod(c(x),c(B))
    }) ~ -1 + eta)$coefficients
    new_llik <- ftr_log_likelihood(input,B,gam_new)
  }
  return(list(gam = gam_old,betas = beta_old,B = compose_parafac(beta_old),llik = llik, total_time = proc.time()[3] - start_model))
}
