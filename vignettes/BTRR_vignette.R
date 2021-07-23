## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Data Generation----------------------------------------------------------
library(bayestensorreg)
set.seed(95064)
simulated_data <-
  TRR_simulated_data(
    subjects = 1,
    max_time = 200,
    margin_sizes = c(50, 50),
    CNR = 1,
    k = 0.3,
    num_active_regions = 2,
    obs.var = 1
  )

## ----True Tensor Coefficient--------------------------------------------------
par(mar=c(0,0,0,0))
image(simulated_data$true_betas, xaxt="n", yaxt="n")

## ---- eval=F------------------------------------------------------------------
#  results <- sapply(seq(3), function(rank) {
#    rank_results <- BTRR_single_subject(
#      input = simulated_data,
#      n.iter = 1000,
#      n.burn = 0,
#      ranks = rank,
#      hyperparameters = NULL,
#      save_llik = T
#    )
#  }, simplify = F)
#  

