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
#' @importFrom parallel detectCores makeCluster parApply stopCluster
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- FTRTucker(input)
#' }
FTRTucker <- function(input,ranks = NULL,epsilon = 1e-4,betas_LASSO = TRUE,num_threads = NULL) {
  start_time <- proc.time()[3]
  Dim <- length(dim(input$X)) - 1
  if(is.null(ranks)) ranks <- rep(1,Dim)
  gam_lm <- lm(input$y ~ -1 + input$eta)
  cat("Log-likelihood without tensor:",logLik(gam_lm),"\n")
  gam_new <- gam_lm$coefficients
  ytil <- c(input$y - input$eta %*% gam_new)
  avail_threads <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(avail_threads)
  B_init <- parallel::parApply(cl,input$X, seq(Dim), function(x, ytil) {
    return(lm(ytil ~ -1 + x)$coefficients)
  }, ytil = ytil)
  parallel::stopCluster(cl)
  beta_new <- sapply(seq(Dim), function(d) {
    b_d <- svd(kFold(B_init,d), nu = ranks[d], nv = ranks[d])$u
    return(b_d)
  }, simplify = F)
  G_new <- sapply(ranks, function(r) seq(r), simplify = F) |>
    expand.grid() |>
    apply(1,function(x) length(unique(x)) == 1) |>
    as.numeric() |>
    array(dim = ranks)
  # beta_new <- mapply(function(p_j,r_j) {
  #   out <- matrix(rnorm(p_j*r_j,sd = 0.025),p_j,r_j)
  #   out[seq(r_j),] <- 1
  #   return(out)
  # },p_j = head(dim(input$X),-1),r_j = ranks,SIMPLIFY = FALSE)
  # G_new <- array(rnorm(prod(ranks)),dim = ranks)
  G_old <- G_new
  beta_old <- beta_new
  gam_old <- gam_new
  llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  new_llik <- max(1000, 5*epsilon)
  step <- 1
  while(abs(new_llik - llik) > epsilon) {
    cat("Step",step,"Log-likelihood",new_llik,"\n")
    G_old <- G_new
    beta_old <- beta_new
    gam_old <- gam_new
    step <- step + 1
    llik <- new_llik
    step_y <- input$y - c(tcrossprod(gam_new,input$eta))
    for(d in seq(length(dim(input$X)) - 1)) {
      # step_X <- t(apply(input$X,length(dim(input$X)),function(X_i){
      #   X_id <- t(apply(X_i,d,identity))
      #   XB_not_d <- X_id %*% Reduce(`%x%`,rev(beta_new[-d])) %*% apply(G_new,d,identity)
      #   return(XB_not_d)
      # }))
      step_X <- kFold(input$X,c(d,length(dim(input$X))))
      step_B <- Reduce(`%x%`,rev(beta_new[-d])) %*% apply(G_new,d,identity)
      # cat("d =",d,"min(B) =",min(c(step_B)),"max(B) =",max(c(step_B)),"\n")
      if(abs(min(c(step_B))) < .Machine$double.eps &
         abs(max(c(step_B))) < .Machine$double.eps) {
        message("Tensor coefficient computationally zero.")
        beta_old[[d]] <- matrix(0,nrow = dim(input$X)[d], ncol = ranks[d])
        return(list(gam = gam_old,B = array(compose_tucker_ftr_vec(beta_old,G_old),
                                            dim = head(dim(input$X),-1)),
                    betas = beta_old, G = G_old, llik = llik,
                    total_time = proc.time()[3] - start_time))
      }
      step_XB <- step_X %*% step_B
      step_X <- matrix(step_XB, nrow = tail(dim(input$X),1), byrow = T)
      # step_X <- kFold(input$X,c(d,length(dim(input$X)))) %*%
      #   Reduce(`%x%`,rev(beta_new[-d])) %*% apply(G_new,d,identity) |>
      #   matrix(nrow = tail(dim(input$X),1), byrow = T)
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
    step_X <- kFold(input$X, length(dim(input$X))) %*% Reduce(`%x%`,rev(beta_new))
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
    BX <- c(kFold(input$X,length(dim(input$X))) %*% compose_tucker_ftr_vec(beta_new,G_new))
    gam_new <- lm(input$y - BX ~ -1 + input$eta)$coefficients
    new_llik <- ftr_log_likelihood(input,compose_tucker_ftr_vec(beta_new,G_new),gam_new)
  }
  return(list(gam = gam_old,B = array(compose_tucker_ftr_vec(beta_old,G_old),
                                      dim = head(dim(input$X),-1)),
              betas = beta_old, G = G_old, llik = llik,
              total_time = proc.time()[3] - start_time))
}
