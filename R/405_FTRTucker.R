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
#'   Tucker tensor decomposition? Defaults to \code{TRUE}.
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
FTRTucker <- function(input,ranks = NULL,epsilon = 1e-3,betas_LASSO = TRUE) {
  start_time <- proc.time()[3]
  if(is.null(ranks)) ranks <- rep(1,length(dim(input$X)) - 1)
  gam_lm <- lm(input$y ~ -1 + input$eta)
  cat("Log-likelihood without tensor:",logLik(gam_lm),"\n")
  gam_new <- gam_lm$coefficients
  beta_new <- mapply(function(p_j,r_j) {
    out <- matrix(rnorm(p_j*r_j,sd = 0.025),p_j,r_j)
    out[seq(r_j),] <- 1
    return(out)
  },p_j = head(dim(input$X),-1),r_j = ranks,SIMPLIFY = FALSE)
  G_new <- array(rnorm(prod(ranks)),dim = ranks)
  ##### Implementing SQUAREM ----
  theta <- c(c(G_new),unlist(beta_new),gam_new)
  # theta_new <- FTRTucker_fixpt(theta, input, ranks, betas_LASSO)
  em_out <- SQUAREM::squarem(
    par = theta,
    fixptfn = FTRTucker_fixpt,
    input = input,
    ranks = ranks,
    betas_LASSO = betas_LASSO,
    control = list(tol = epsilon ^ (1 / length(theta)), trace = TRUE)
  )
  theta_out <- untheta(em_out$par,ranks = ranks, dim_X = head(dim(input$X),-1))
  ##### End Implementing SQUAREM ----
  G_old <- G_new
  beta_old <- beta_new
  gam_old <- gam_new
  llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  new_llik <- llik + 5*epsilon
  step <- 1
  while(new_llik - llik > epsilon) {
    cat("Step",step,"Log-likelihood",new_llik,"\n")
    G_old <- G_new
    beta_old <- beta_new
    gam_old <- gam_new
    step <- step + 1
    llik <- new_llik
    step_y <- input$y - c(tcrossprod(gam_new,input$eta))
    for(d in seq(length(dim(input$X)) - 1)) {
      step_X <- t(apply(input$X,length(dim(input$X)),function(X_i){
        X_id <- t(apply(X_i,d,identity))
        XB_not_d <- X_id %*% Reduce(`%x%`,rev(beta_new[-d])) %*% apply(G_new,d,identity)
        return(XB_not_d)
      }))
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
    step_X <- t(apply(input$X,length(dim(input$X)),function(X_i) t(Reduce(`%x%`,rev(beta_new)))%*%c(X_i)))
    if(nrow(step_X) == 1) step_X <- t(step_X)
    if(length(G_new) < 2){
      G_new <- array(lm(step_y ~ -1 + step_X)$coefficients,dim = ranks)
    }else{
      cv.G_new <- try(glmnet::cv.glmnet(step_X,step_y,type.measure = "mse",alpha = 1,
                            family = "gaussian",intercept = FALSE))
      if(class(cv.G_new) == "cv.glmnet") {
          G_new <- array(c(as.matrix(coef(cv.G_new,s = "lambda.min"))[-1,]),
                     dim = ranks)
      }else{
        G_new <- array(c(lm(step_y ~ -1 + step_X)$coefficients), dim = ranks)
      }
    }
    gam_new <- lm(input$y - apply(input$X,length(dim(input$X)),function(X_i){
      crossprod(c(X_i),c(compose_tucker_ftr_vec(beta_new,G_new)))
    }) ~ -1 + input$eta)$coefficients
    new_llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  }
  return(list(gam = gam_old,B = array(compose_tucker_ftr_vec(beta_old,G_old),
                                      dim = head(dim(input$X),-1)),
              betas = beta_old, G = G_old, llik = llik,
              total_time = proc.time()[3] - start_time))
}

untheta <- function(theta, ranks, dim_X) {
  G_inds <- seq(prod(ranks))
  beta_inds <- seq(sum(dim_X * ranks)) + max(G_inds)
  gam_inds <- setdiff(seq(length(theta)), c(G_inds,beta_inds))
  G_out <- array(theta[G_inds], dim = ranks)
  theta_remain <- theta[-G_inds]
  beta_out <- vector("list", length = length(ranks))
  for(j in 1:length(ranks)) {
    beta_out[[j]] <-
      matrix(theta_remain[seq(dim_X[j] * ranks[j])], dim_X[j], ranks[j])
    theta_remain <- theta_remain[-seq(dim_X[j] * ranks[j])]
  }
  return(
    list(G = G_out, beta = beta_out, gam = theta_remain)
  )
}

FTRTucker_fixpt <- function(theta, input, ranks, betas_LASSO) {
  model_params <- untheta(theta, ranks, head(dim(input$X),-1))
  step_y <- input$y - c(tcrossprod(model_params$gam, input$eta))
  for (d in seq(length(dim(input$X)) - 1)) {
    step_X <- t(apply(input$X, length(dim(input$X)), function(X_i) {
      X_id <- t(apply(X_i, d, identity))
      XB_not_d <-
        X_id %*% Reduce(`%x%`, rev(model_params$beta[-d])) %*% apply(model_params$G, d, identity)
      return(XB_not_d)
    }))
    if (betas_LASSO) {
      cv.beta_new <-
        try(glmnet::cv.glmnet(
          step_X,
          step_y,
          type.measure = "mse",
          alpha = 1,
          family = "gaussian",
          intercept = FALSE
        ))
      if (class(cv.beta_new) == "cv.glmnet") {
        model_params$beta[[d]] <-
          matrix(c(as.matrix(coef(
            cv.beta_new, s = "lambda.min"
          ))[-1, ]),
          ncol = ranks[d])
      } else{
        model_params$beta[[d]] <- matrix(lm(step_y ~ -1 + step_X)$coefficients,
                                ncol = ranks[d])
      }
    } else{
      model_params$beta[[d]] <- matrix(lm(step_y ~ -1 + step_X)$coefficients,
                              ncol = ranks[d])
    }
    model_params$beta[[d]][is.na(model_params$beta[[d]])] <-
      0 # In some cases with MRI data, all scan values are equal to zero.
    model_params$beta[[d]][seq(ranks[d]), ] <- 1
  }
  step_X <-
    t(apply(input$X, length(dim(input$X)), function(X_i)
      t(Reduce(
        `%x%`, rev(model_params$beta)
      )) %*% c(X_i)))
  if (nrow(step_X) == 1)
    step_X <- t(step_X)
  if (length(G_new) < 2) {
    G_new <- array(lm(step_y ~ -1 + step_X)$coefficients, dim = ranks)
  } else{
    cv.G_new <-
      try(glmnet::cv.glmnet(
        step_X,
        step_y,
        type.measure = "mse",
        alpha = 1,
        family = "gaussian",
        intercept = FALSE
      ))
    if (class(cv.G_new) == "cv.glmnet") {
      model_params$G <-
        array(c(as.matrix(coef(cv.G_new, s = "lambda.min"))[-1,]),
              dim = ranks)
    } else{
      model_params$G <-
        array(c(lm(step_y ~ -1 + step_X)$coefficients), dim = ranks)
    }
  }
  model_params$gam <-
    lm(input$y - apply(input$X, length(dim(input$X)), function(X_i) {
      crossprod(c(X_i), c(compose_tucker_ftr_vec(model_params$beta, model_params$G)))
    }) ~ -1 + input$eta)$coefficients
  return(unlist(model_params))
}
