#' Frequentist tensor regression with the Tucker decomposition
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param ranks The ranks to be used with the Tucker decomposition. This should
#'   be a vector with the same length as the tensor covariate for each subject.
#' @param epsilon a value for the stopping rule of the algorithm. Specifically,
#'   this is the upper bound for the differences in the log-likelihood between
#'   two iterations of the algorithm.
#' @param betas_LASSO (logical) Should the LASSO be applied to the betas in the
#'   Tucker tensor decomposition? Defaults to \code{FALSE}.
#' @param G_LASSO (logical) Should the LASSO be applied to the core tensor in the
#'   Tucker tensor decomposition? Defaults to \code{TRUE}.
#' @param step_limit The maximum number of steps that can be taken before
#'   deciding that the algorithm did not converge
#'
#' @importFrom stats lm rnorm logLik coef
#' @import glmnet
#'
#' @return  A list with elements \code{gam} (vector coefficient result),
#'   \code{betas} (tensor decomposition components), \code{G} (the core tensor in the tensor decomposition),
#'   \code{B} (the tensor coefficient), \code{llik} (the value of the
#'   log-likelihood) and \code{total_time} (time spent to complete the analysis).
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- FTRTucker(input)
#' }
FTRTucker <-
  function(
    input,
    ranks = NULL,
    epsilon = 1e-4,
    betas_LASSO = FALSE,
    G_LASSO = TRUE,
    step_limit = 1000
    ) {
  start_time <- proc.time()[3]
  Dim <- length(dim(input$X))-1
  sample_size <- length(input$y)
  if(is.null(ranks)) ranks <- rep(1,Dim)
  if(!is.null(input$eta)) {
    gam_lm <- lm(input$y ~ -1 + input$eta)
    llik <- c(logLik(gam_lm))
    gam_new <- gam_lm$coefficients
  } else {
    gam_new <- NULL
    llik <- sum(dnorm(input$y, mean = mean(input$y), sd = sd(input$y), log = TRUE))
  }
  cat("Log-likelihood without tensor:",llik,"\n")

  beta_new <- mapply(function(p_j,r_j) {
    out <- matrix(rnorm(p_j*r_j,sd = 0.025),p_j,r_j)
    out[seq(r_j),] <- 1
    return(out)
  },p_j = head(dim(input$X),-1),r_j = ranks,SIMPLIFY = FALSE)
  G_new <- array(rnorm(prod(ranks)),dim = ranks)
  new_llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  step <- 1
  while(100*abs(new_llik - llik)/abs(llik) > epsilon & step < step_limit) {
    cat("Step",step,"Log-likelihood",new_llik,"\n")
    G_old <- G_new
    beta_old <- beta_new
    gam_old <- gam_new
    step <- step + 1
    llik <- new_llik
    if(!is.null(input$eta)) {
      step_y <- input$y - c(tcrossprod(gam_new,input$eta))
    } else {step_y <- input$y}
    for(d in seq(Dim)) {
      # step_X_old <- t(apply(input$X,length(dim(input$X)),function(X_i){
      #   X_id <- t(apply(X_i,d,identity))
      #   XB_not_d <- X_id %*% Reduce(`%x%`,rev(beta_new[-d])) %*% apply(G_new,d,identity)
      #   return(XB_not_d)
      # }))
      step_B <- tcrossprod(Reduce(`%x%`,rev(beta_new[-d])),kFold(G_new,d))
      step_X <- kFold(input$X,c(Dim+1,d)) %*% step_B |>
        matrix(nrow = sample_size, byrow = F)
      if(betas_LASSO){
        cv.beta_new <- try(glmnet::cv.glmnet(step_X,step_y,type.measure = "mse",alpha = 1,
                                 family = "gaussian",intercept = FALSE))
        if(class(cv.beta_new) == "cv.glmnet") {
          beta_new[[d]] <- matrix(c(as.matrix(coef(cv.beta_new,s = "lambda.min"))[-1,]),
                                  ncol = ranks[d])
        }else{
          beta_new[[d]] <- matrix(lm(step_y ~ -1 + step_X)$coefficients,
                                  ncol = ranks[d])
        }
        # beta_new[[d]] <- matrix(c(as.matrix(coef(cv.beta_new,s = "lambda.min"))[-1,]),ncol = ranks[d])
      }else{
        beta_new[[d]] <- matrix(lm(step_y ~ -1 + step_X)$coefficients,
                                ncol = ranks[d])
      }
      beta_new[[d]][is.na(beta_new[[d]])] <- 0 # In some cases with MRI data, all scan values are equal to zero.
      beta_new[[d]][seq(ranks[d]),] <- 1
    }
    # step_X <- t(apply(input$X,length(dim(input$X)),function(X_i) t(Reduce(`%x%`,rev(beta_new)))%*%c(X_i)))
    step_X <- kFold(input$X, Dim + 1) %*% Reduce(`%x%`,rev(beta_new))
    if(nrow(step_X) == 1) step_X <- t(step_X)
    if(length(G_new) < 2){
      G_new <- array(lm(step_y ~ -1 + step_X)$coefficients,dim = ranks)
    }else{
      if(G_LASSO){
        cv.G_new <- try(glmnet::cv.glmnet(step_X,step_y,type.measure = "mse",alpha = 1,
                                          family = "gaussian",intercept = FALSE))
        if(class(cv.G_new) == "cv.glmnet") {
          G_new <- array(c(as.matrix(coef(cv.G_new,s = "lambda.min"))[-1,]),
                         dim = ranks)
        }else{
          G_new <- array(c(lm(step_y ~ -1 + step_X)$coefficients), dim = ranks)
        }
      }else{
        G_new <- array(c(lm(step_y ~ -1 + step_X)$coefficients), dim = ranks)
      }
    }
    BX <- c(kFold(input$X,Dim + 1) %*% compose_tucker_ftr_vec(beta_new,G_new))
    if(!is.null(input$eta)) {
      gam_new <- lm(input$y - BX ~ -1 + input$eta)$coefficients
    }
    new_llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  }
  return(list(gam = gam_old,B = array(compose_tucker_ftr_vec(beta_old,G_old),
                                      dim = head(dim(input$X),-1)),
              betas = beta_old, G = G_old, llik = llik,
              total_time = proc.time()[3] - start_time))
}
