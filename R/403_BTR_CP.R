#' Bayesian tensor regression with the CP decomposition
#'
#' This code applies the model by Guhaniyogi et. al. (2017)
#'
#' @param input An object of class \code{TR_data} that contains (at least) the
#'   elements \code{y} (a vector of response values) and \code{X} (an array of
#'   covariate values). Optionally, \code{eta} (a matrix of nuisance covariates)
#'   can also be included. Other list elements will be ignored.
#' @param max_rank (a scalar) The rank of the CP decomposition to be used
#' @param n_iter  (a scalar) the number of posterior samples desired
#' @param n_burn (a scalar) the number of posterior samples to discard as a
#'   burn-in
#' @param hyperparameters  a list with the (scalar) elements \code{a.tau},
#'   \code{b.tau}, \code{a.lambda}, \code{b.lambda}, \code{a.epsilon}, and
#'   \code{b.epsilon} defining the values of the hyperparameters within the
#'   model. If \code{NULL}, then default values will be used.
#' @param save_dir (a character) A path to a directory in which the temporary
#'   results will be saved. Defaults to \code{NULL}. If \code{NULL}, then
#'   temporary saves are not made.
#' @param num_threads the number of threads that can be used to find initial
#'   conditions
#'
#' @importFrom stats rnorm dgamma
#'
#' @return A list with the posterior samples
#' @export
#'
#' @examples
#' \dontrun{
#' input <- TR_simulated_data()
#' results <- BTR_CP(input)
#' }
BTR_CP <- function(input,max_rank = 1,n_iter = 100, n_burn = 0,hyperparameters = NULL, save_dir = NULL, num_threads = NULL) {
  # Logic checks for data compatibility
  if(length(input$eta) == 0) input$eta <- matrix(0,nrow = length(input$y), ncol = 1)
  if(sum(as.numeric(c("y","X","eta") %in% names(input))) != 3) stop("Input must be a list with at least the elements y, X, and eta.")

  # Summary values pulled from the data
  p <- head(dim(input$X),-1)
  D <- length(p)

  # Set hyperparameters ----
  if(is.null(hyperparameters)) {
    a_lam <- 3 # From Guhaniyogi et al. [2017]
    b_lam <- a_lam^(1/(2*D)) # From Guhaniyogi et al. [2017]
    nu <- 2 # From Guhaniyogi et al. [2017]
    s_02 <- -log(0.95)
    a_sig <- 3
    b_sig <- 20
    invSig_0 <- (1/900)*diag(ncol(input$eta)) # From Guhaniyogi et al. [2017]
    mu_gam = rep(0,ncol(input$eta)) # From Guhaniyogi et al. [2017]
    alpha_grid <- seq(max_rank^(-D),
                      max_rank^(-.1),
                      length.out = 10) # From Guhaniyogi et al. [2017]
    a_tau <- 1 # This is what was in Shaan's code
    b_tau <- max_rank^(1/D - 1) # From Guhaniyogi et al. [2017]
  }else{
    if(!is.list(hyperparameters)) stop("hyperparameters should be NULL or a list with the elements a_lam (a scalar), b_lam (a scalar), a_sig (a scalar), b_sig (a scalar), invSig_0 (a square matrix with dimension equal to the number of columns in eta), alpha_grid (a vector), a_tau (a scalar), and b_tau (a scalar)")
    if(sum(as.numeric(c("a_lam","b_lam","nu","s_02","invSig_0","mu_gam","alpha_grid","a_tau","b_tau") %in% names(hyperparameters))) != 9) stop("hyperparameters should be NULL or a list with the elements a_lam (a scalar), b_lam (a scalar), a_sig (a scalar), b_sig (a scalar), invSig_0 (a square matrix with dimension equal to the number of columns in eta), alpha_grid (a vector), a_tau (a scalar), and b_tau (a scalar)")
    list2env(hyperparameters,globalenv())
  }

  # Set storage ----
  results <-
    list(
      betas = list(),
      W = list(),
      tau = numeric(n_iter),
      lambda = list(),
      Phi = matrix(NA,n_iter,max_rank),
      gam = matrix(NA,n_iter,ncol(input$eta)),
      sig_y2 = numeric(n_iter),
      llik = numeric(n_iter)
    )

  # Set initials ----
  # betas <- sapply(seq(D),function(j){
  #   matrix(rnorm(p[j]*max_rank,sd = 0.025),p[j],max_rank)
  # },simplify = FALSE)
  gam <- lm(input$y ~ -1 + input$eta)$coefficients
  if(any(is.na(gam))) gam <- 0
  ytil <- c(input$y - input$eta %*% gam)
  avail_threads <- parallel::detectCores() - 1
  if(is.null(num_threads)) num_threads <- avail_threads
  num_threads <- min(num_threads, avail_threads)
  cl <- parallel::makeCluster(num_threads)
  B_init <- parallel::parApply(cl,input$X, seq(D), function(x, ytil) {
    return(lm(ytil ~ -1 + x)$coefficients)
  }, ytil = ytil)
  parallel::stopCluster(cl)
  betas <- sapply(seq(D), function(d) {
    b_d <- svd(kFold(B_init,d), nu = max_rank, nv = max_rank)$u
    return(b_d)
  }, simplify = F)
  lam <- sapply(seq(D),function(d) rep(1,max_rank), simplify = FALSE)
  tau <- 1
  sig_y2 <- 1
  W <- sapply(p,function(pp) matrix(1,pp,max_rank), simplify = FALSE)
  M <- 4 # The number of samples for each alpha_grid value in the Griddy Gibbs

  # Begin MCMC ----
  start_time <- proc.time()[3]
  for(s in seq(n_iter)){
    # Griddy Gibbs to draw a value for alpha ----
    if(max_rank > 1){
      Cr <- t(sapply(seq(D),function(j){sapply(seq(max_rank),function(r){crossprod(betas[[j]][,r],diag(1/W[[j]][,r]))%*%betas[[j]][,r]})})) # D x R
      weights <- sapply(alpha_grid,function(alpha){
        Phi_l <- cp_draw_Phi(M,alpha,betas,W,b_tau) # M x R
        tau_l <- cp_draw_tau(M,a_tau,b_tau,betas,W,Phi_l) # length M
        ldPhi <- apply(Phi_l,1,function(phi_l) {
          (lgamma(alpha*max_rank) - sum(lgamma(rep(alpha,max_rank)))) +
            sum((alpha - 1)*log(phi_l))
        })
        ldtau <- dgamma(tau_l,alpha*max_rank,b_tau,log = T)
        tau_r <- apply(Phi_l,2,`*`,tau_l)
        ldbetas <- apply(tau_r,1,function(tr){sum(-(sum(p)/2)*log(tr) - colSums(Cr)/tr)})
        out <- ldPhi + ldtau + ldbetas
      })
      max_weight <- max(weights)
      alpha_probs <- apply(weights,2,function(w) mean(exp(w - max_weight)))
      alpha <- sample(alpha_grid,size=1,prob = alpha_probs)
    }else{
      alpha <- 1
    }

    # Draw Phi ----
    Phi <- cp_draw_Phi(1,alpha,betas,W,b_tau)
    # Draw tau ----
    tau <- cp_draw_tau(1,max_rank*alpha,alpha*(max_rank/nu)^(1/D),betas,W,Phi)
    # Draw lambda ----
    lam <- cp_draw_lam(a_lam,b_lam,betas,Phi,tau)
    # Draw W ----
    W <- cp_draw_W(lam,betas,Phi,tau)
    # Draw betas ----
    for(r in seq(max_rank)) {
      for(j in seq(D)) {
        betas[[j]][,r] <- cp_draw_betas(X = input$X,betas = betas,
                                        sig_y2 = sig_y2,y = input$y,
                                        gam = gam,eta = input$eta,
                                        W = W,Phi = Phi,tau = tau,j = j,r = r)
      }
    }
    # Draw sig_y2 ----
    B_s <- btr_compose_parafac(betas)
    XB_s <- c(kFold(input$X,D + 1) %*% c(B_s))
    y_til <- input$y -
      # apply(input$X,length(dim(input$X)),function(X_i){
      #   crossprod(c(X_i),c(B_s))
      # }) -
      XB_s -
      c(input$eta %*% gam)
    sig_y2 <- cp_draw_sig_y2_two(a_sig,b_sig, y_til)
    # Draw gam ----
    if(all(c(input$eta) != 0)){
      y_til <- input$y - XB_s
        # apply(input$X,length(dim(input$X)),function(X_i){
        #   crossprod(c(X_i),c(B_s))
        # })
      gam <- cp_draw_gam_two(input$eta,invSig_0,mu_gam,y_til,sig_y2)
    }
    # Calculate log-likelihood ----
    llik <- sum(dnorm(
      input$y - XB_s -
        # apply(input$X,
        #       length(dim(input$X)),
        #       function(X_i)
        #         crossprod(c(X_i), c(B_s))) -
        input$eta %*% gam,
      mean = 0,
      sd = sqrt(sig_y2),
      log = TRUE
    ))

    # Store results ----
    results$betas[[s]] <- betas
    results$W[[s]] <- W
    results$tau[s] <- tau
    results$lambda[[s]] <- lam
    results$Phi[s,] <- Phi
    results$gam[s,] <- gam
    results$sig_y2[s] <- sig_y2
    results$llik[s] <- llik

    # Status report ----
    if(s %% floor(n_iter / 10) == 0) cat("Iteration",s,"of",n_iter,"Log-likelihood:",llik,"\n")
    if(!is.null(save_dir) & s %% floor(n_iter / 10) == 0) {
      results$total_time <- proc.time()[3] - start_time
      saveRDS(results,file.path(save_dir,paste0("BTR_cp_y_X_temp_rank_",paste0(max_rank,collapse = ""),"_",format(Sys.Date(),"%Y%m%d"),".rds")))
    }
  }
  # Discard burn-in ----
  if(n_burn > 0){
    results$betas <- results$betas[-seq(n_burn)]
    results$W <- results$W[-seq(n_burn)]
    results$tau <- results$tau[-seq(n_burn)]
    results$lambda <- results$lambda[-seq(n_burn)]
    results$Phi <- results$Phi[-seq(n_burn),]
    results$gam <- results$gam[-seq(n_burn),]
    results$sig_y2 <- results$sig_y2[-seq(n_burn)]
  }
  results$total_time <- proc.time()[3] - start_time
  return(results)
}
