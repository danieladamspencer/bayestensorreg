#' Simulated tensor regression data
#'
#' @param subjects A positive integer indicating the number of subjects for
#'   which to generate data
#' @param tensor_dims A vector that indicates the dimensions of the tensor
#'   covariate for each subject
#' @param CNR The contrast-to-noise ratio, defined by Welvaert and Rosseel
#'   (2013) as the size of the (largest) true coefficient value divided by
#'   the observational noise.
#' @param num_active The number of contiguous nonzero activation regions
#' @param other_covar A vector of true values for other (nuisance) coefficients
#'   that will be used to generate non-tensor-valued covariates
#'
#' @return An object of class \code{TR_data} that contains the
#'   elements \code{y} (a vector of response values), \code{X} (an array of
#'   covariate values), \code{eta} (a matrix of nuisance covariates),
#'   \code{true_B} (an array with the true tensor coefficient values), and
#'   \code{gam} (a vector with the true nuisance coefficient values)
#' @export
#' @importFrom neuRosim specifyregion
#'
#' @examples
#' input <- TR_simulated_data()
TR_simulated_data <-
  function(subjects = 400,
           tensor_dims = c(50, 50),
           CNR = 1,
           num_active = 3,
           other_covar = c(1, 1)) {
    B <- Reduce(`+`, sapply(seq(num_active), function(zz) {
      neuRosim::specifyregion(
        dim = tensor_dims,
        coord = sapply(tensor_dims, function(z)
          sample(seq(z), size = 1)),
        radius = ceiling(min(tensor_dims) * runif(1, 0.02, 0.1)),
        form = "sphere",
        fading = runif(1, 0.5, 1)
      )
    }, simplify = FALSE))
    eta <-
      matrix(rnorm(subjects * length(other_covar)), subjects, length(other_covar))
    gam <- other_covar
    X <-
      array(rnorm(prod(tensor_dims) * subjects), dim = c(tensor_dims, subjects))
    y <-
      apply(X, length(dim(X)), function(xx)
        sum(xx * B * CNR)) + c(eta %*% gam) + rnorm(subjects)
    return(list(
      y = y,
      X = X,
      true_B = B,
      eta = eta,
      gam = gam
    ))
  }

#' Simulated data for tensor response regression with Gaussian graphical model
#'
#' @param subjects The number of subjects in the data (a scalar)
#' @param regions The number of response tensors per subject in the simulated
#'   data
#' @param max_time The length of the tensor time series (a scalar)
#' @param margin_sizes The sizes of the tensor coefficients (a vector)
#' @param SNR The signal-to-noise ratio (a scalar), as defined in
#'   Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise
#'   ratio and contrast-to-noise ratio for fMRI data. PloS one, 8(11), e77089.
#' @param CNR The contrast-to-noise ratio (a scalar), as defined in
#'   Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise
#'   ratio and contrast-to-noise ratio for fMRI data. PloS one, 8(11), e77089.
#' @param conn_regions The number of response tensors that should have nonzero
#'   partial correlations
#' @param conn_level The partial correlation for the connected regions
#'
#' @return A list with class \code{TRR_GGM_data} with elements \code{Y},
#'   \code{x}, \code{true_d}, \code{true_B}, and \code{true_d_covar}
#' @export
#'
#' @importFrom neuRosim specifyregion specifydesign
#'
#' @examples
#' \dontrun{
#' input <- TRR_GGM_simulated_data()
#' }
TRR_GGM_simulated_data <-
  function(subjects = 20,
           regions = 10,
           max_time = 200,
           margin_sizes = c(10, 10, 10),
           SNR = 1,
           CNR = 0.05,
           conn_regions = 2,
           conn_level = 0.9) {
    # Construct covariance based on connectivity
    d_covar <- diag(regions)
    conns <- matrix(NA, conn_regions, 2)
    for (g in seq(conn_regions)) {
      conns[g, ] <- sample(seq(regions), size = 2, replace = FALSE)
      if (g >= 2) {
        for (gg in seq(g - 1)) {
          while (identical(conns[g, ], conns[gg, ]) |
                 identical(conns[g, ], rev(conns[gg, ]))) {
            conns[g, ] <- sample(seq(regions), size = 2, replace = FALSE)
          }
        }
      }
      d_covar[conns[g, 1], conns[g, 2]] <-
        d_covar[conns[g, 2], conns[g, 1]] <- conn_level
    }
    # Create the subject-ROI effects
    d <-
      matrix(rnorm(subjects * regions), subjects, regions) %*% chol(d_covar) * sqrt(SNR)
    # Make the coefficient tensors
    p <- sapply(seq(regions), function(g)
      margin_sizes, simplify = F)
    B <- sapply(p, function(region_idx) {
      neuRosim::specifyregion(
        dim = region_idx,
        coord = runif(length(region_idx)) * region_idx,
        radius = min(round(min(region_idx * 0.1)),
                     sample(3:5, 1)),
        form = "sphere"
      ) * CNR
    }, simplify = FALSE)

    # Create a task covariate
    x <-
      neuRosim::specifydesign(
        onsets = seq(0.1 * max_time, 0.9 * max_time, length.out = 5),
        durations = 1,
        totaltime = max_time,
        TR = 1,
        effectsize = 1,
        conv = "double-gamma",
        param = list(list(
          a1 = 6,# Delay of response relative to onset
          a2 = 12,# Delay of undershoot relative to onset
          b1 = 0.9,# Dispersion of response
          b2 = 0.9,# Dispersion of undershoot
          c = 0.35 # Scale of undershoot
        ))
      )
    X <- matrix(x, nrow = max_time, ncol = subjects)
    # Make the response data
    Y <- mapply(function(b, dd) {
      b %o% X +
        sapply(dd,
               array,
               dim = c(dim(b), max_time),
               simplify = "array") +
        array(rnorm(prod(dim(b)) * max_time * subjects),
              dim = c(dim(b), max_time, subjects))
    },
    b = B,
    dd = split(d, col(d)),
    SIMPLIFY = FALSE)

    output <-
      list(
        Y = Y,
        x = X,
        true_d = d,
        true_B = B,
        true_d_covar = d_covar
      )

    class(output) <- "TRR_GGM_data"
    return(output)
  }

#' Simulated tensor response regression data
#'
#' @param subjects The number of subjects in the data (a scalar)
#' @param max_time The length of the tensor time series (a scalar)
#' @param margin_sizes The sizes of the tensor coefficient (a vector)
#' @param CNR The contrast-to-noise ratio (a scalar), as defined in
#'   Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise
#'   ratio and contrast-to-noise ratio for fMRI data. PloS one, 8(11), e77089.
#' @param k The AR(1) autoregression coefficient (a scalar)
#' @param num_active_regions The number of nonzero regions in the tensor coefficient
#' @param obs.var The observation variance (a scalar)
#'
#' @importFrom neuRosim specifyregion specifydesign
#'
#' @return A list with objects \code{Y} (the tensor-valued response), \code{x} (the matrix-valued covariate), \code{true_betas} (the array of true tensor coefficient values), and \code{true_k} the true value for the autocorrelation coefficient.
#' @export
#'
#' @examples
#' input <- TRR_simulated_data()
TRR_simulated_data <-
  function(subjects = 1,
           max_time = 200,
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
        onsets = seq(0.1 * max_time, 0.9 * max_time, length.out = 5),
        durations = 1,
        totaltime = max_time,
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
    Sigma_AR <- toeplitz(k ^ seq(0, max_time - 1)) * obs.var / (1 - k ^ 2)

    # Find the expected value of the response at all points
    true_mean <- betas %o% c(x)

    # Simulate the response tensor(s)
    y <- sapply(seq(subjects), function(each_subject) {
      ar_series <- apply(betas, seq(length(dim(betas))), function(b) {
        error_series <- chol(Sigma_AR) %*% rnorm(max_time, sd = sqrt(obs.var))
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
      x = matrix(x, max_time, subjects),
      true_betas = betas,
      k = k
    )

    # Assign a class to ensure compatibility
    class(simulated_data) <- "TRR_data"

    return(simulated_data)
  }

