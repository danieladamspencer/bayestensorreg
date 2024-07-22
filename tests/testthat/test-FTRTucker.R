test_that("FTRTucker works without eta", {
  set.seed(123)
  input <- TR_simulated_data(subjects = 10, tensor_dims = c(20, 20), other_covar = NULL)
  ranks = c(3,3)
  epsilon = 1e-4
  betas_LASSO = FALSE
  G_LASSO = TRUE
  step_limit = 1000
  model_output <- FTRTucker(input = input, ranks = c(1,1))
  expect_true(is.null(model_output$gam))
})
