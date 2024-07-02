test_that("FTRTucker works without eta", {
  set.seed(123)
  input <- TR_simulated_data(subjects = 20, tensor_dims = c(20, 20), other_covar = NULL)
  model_output <- FTRTucker(input = input, ranks = c(1,1))
  expect_true(is.null(model_output$gam))
})
