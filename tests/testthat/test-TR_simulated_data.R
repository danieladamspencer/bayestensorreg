test_that("simulation works without eta", {
  set.seed(123)
  sim_data <- TR_simulated_data(subjects = 20, tensor_dims = c(20, 20), CNR = 1, num_active = 1, other_covar = NULL)
  expect_true(is.null(sim_data$eta))
  expect_true(is.null(sim_data$gam))
})
