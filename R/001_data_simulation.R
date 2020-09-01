#' Simulated tensor response regression data for a single subject
#'
#' @param subjects The number of subjects in the data (a scalar)
#' @param n.time The length of the tensor time series (a scalar)
#' @param margin_sizes The sizes of the tensor coefficient (a vector)
#' @param CNR The contrast-to-noise ratio (a scalar), as defined in Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise ratio and contrast-to-noise ratio for fMRI data. PloS one, 8(11), e77089.
#' @param k The AR(1) autoregression coefficient (a scalar)
#' @param num_active_regions The number of nonzero regions in the tensor coefficient
#' @param obs.var The observation variance (a scalar)
#'
#' @importFrom neuRosim specifyregion specifydesign
#' @import stats
#'
#' @return A list with objects \code{Y} (the tensor-valued response), \code{x} (the matrix-valued covariate), \code{true_betas} (the array of true tensor coefficient values), and \code{true_k} the true value for the autocorrelation coefficient.
#' @export
#'
#' @examples
#' input <- TRR_simulated_data()
TRR_simulated_data <-
  function(subjects = 1,
           n.time = 200,
           margin_sizes = c(20, 20),
           CNR = 1,
           k = 0.3,
           num_active_regions = 1,
           obs.var = 1) {
    requireNamespace("neuRosim")
    # Create a signal region
    betas <- sapply(seq(num_active_regions), function(nr) {
      active_image <- neuRosim::specifyregion(dim = margin_sizes,
                                    coord = margin_sizes * runif(length(margin_sizes)),
                                    radius = min(round(min(margin_sizes) *
                                                         0.1), sample(3:5, 1)),
                                    form = "sphere", fading = 0.5)
      return(active_image)
    },simplify = F)
    betas <- Reduce(`+`,betas)
    # Scale the active region to correspond to the contrast-to-noise ratio
    betas <- betas * CNR * sqrt(obs.var)
    # Create a task covariate
    x <-
      neuRosim::specifydesign(
        onsets = seq(0.1 * n.time, 0.9 * n.time, length.out = 5),
        durations = 1,
        totaltime = n.time,
        TR = 1,
        effectsize = 1,
        conv = "double-gamma",
        param = list(list(
          a1 = 6,# Delay of response relative to onset
          a2 = 12,# Delay of undershoot relative to onset
          b1 = 0.9,# Dispersion of response
          b2 = 0.9,# Dispersion of undershoot
          c = 0.15 # Scale of undershoot
        ))
      )

    # Create the autoregressive covariance
    Sigma_AR <- toeplitz(k ^ seq(0, n.time - 1)) * obs.var / (1 - k ^ 2)

    # Find the expected value of the response at all points
    true_mean <- betas %o% c(x)

    # Simulate the response tensor(s)
    y <- sapply(seq(subjects), function(each_subject) {
      ar_series <- apply(betas, seq(length(dim(betas))), function(b) {
        error_series <- chol(Sigma_AR) %*% rnorm(n.time, sd = sqrt(obs.var))
        return(error_series)
      })
      correct_dim <-
        c(tail(seq(length(dim(
          ar_series
        ))), -1), head(seq(length(dim(
          ar_series
        ))), 1))
      subject_data <- aperm(ar_series, correct_dim) + true_mean
      return(subject_data)
    }, simplify = "array")

    # Bring all of the data together into a list
    simulated_data <- list(
      Y = y,
      x = matrix(x, n.time, subjects),
      true_betas = betas,
      k = k
    )

    # Assign a class to ensure compatibility
    class(simulated_data) <- "TRR_data"

    return(simulated_data)
  }
