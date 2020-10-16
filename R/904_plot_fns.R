#' Time series ggplot
#'
#' @param x A numeric vector
#' @param title A character string of the plot title
#'
#' @import ggplot2
#'
#' @return A ggplot grob
#' @export
ggplot_ts <- function(x, title = NULL) {
  Time <- NULL
  ts_df <- as.data.frame(cbind(Time = seq(length(x)), x))
  out <- ggplot(ts_df) +
    geom_line(aes(x = Time, y = x)) +
    ggtitle(ifelse(is.null(title), "", title)) +
    theme_bw()
  return(out)
}

#' Plot posterior densities of the tensor coefficient
#'
#' @param B An array in which the last index corresponds to the \eqn{i}th draw
#'   from the posterior distribution.
#' @param truth An optional array of true values for the tensor coefficient
#'   which will be used to separate the density plots for voxels with true zero
#'   values and true nonzero plots.
#'
#' @importFrom reshape2 melt
#' @importFrom tidyr unite
#' @import ggplot2
#'
#' @return A ggplot grob
#' @export
B_densities <- function(B, truth = NULL) {
  voxel <- value <- NULL
  D <- length(dim(B)) - 1
  if (is.null(truth)) {
    coef_df <-
      reshape2::melt(B, varnames = c(paste0("Var", seq(D)), "iter"))
    voxelwise_df <-
      tidyr::unite(coef_df, voxel, grep("Var", names(coef_df), value = TRUE))
    out <- ggplot2::ggplot(voxelwise_df) +
      ggplot2::geom_density(ggplot2::aes(x = value, group = voxel)) +
      ggplot2::labs(x = "", y = "Posterior Density") +
      ggplot2::theme_bw()
  } else {
    if (!identical(dim(truth), head(dim(B), -1)))
      stop("The dimensions of truth should match all but the last dimensions of B.")
    coef_df <-
      reshape2::melt(B, varnames = c(paste0("Var", seq(D)), "iter"))
    truth_df <-
      reshape2::melt(truth, value.name = "Truth")
    truth_df$Truth <-
      factor(truth_df$Truth != 0, labels = paste("True", c("Zero", "Nonzero")))
    coef_truth_df <- merge(coef_df, truth_df)
    voxelwise_df <-
      tidyr::unite(coef_truth_df, voxel, grep("Var", names(coef_truth_df),
                                              value = TRUE))
    out <- ggplot2::ggplot(voxelwise_df) +
      ggplot2::geom_density(ggplot2::aes(x = value, group = voxel)) +
      ggplot2::facet_grid(Truth ~ .) +
      ggplot2::labs(x = "", y = "Posterior Density") +
      ggplot2::theme_bw()
  }
  return(out)
}

#' Make a tile plot using ggplot2
#'
#' @param x a matrix
#' @param title the plot title
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @return a ggplot grob
#' @export
gg_tile_plot <- function(x,title = NULL) {
  if(!is.matrix(x)) stop("The input to this function should be a matrix.")
  Var1 <- Var2 <- value <- NULL
  x_df <- reshape2::melt(x)
  out <- ggplot(x_df) +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    scale_color_gradient2("") +
    scale_fill_gradient2("") +
    ggtitle(ifelse(is.null(title),"",title)) +
    labs(x="",y="") +
    theme_bw()
  return(out)
}
