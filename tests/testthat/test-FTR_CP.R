test_that("Works with no eta", {
  set.seed(123)
  input <- TR_simulated_data(subjects = 20, tensor_dims = c(20, 20), CNR = 1, num_active = 1, other_covar = NULL)
  model_output <- FTR_CP(input = input, rank = 1)
  expect_true(is.null(model_output$gam))
})
