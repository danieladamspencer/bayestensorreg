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

#' Melt a matrix into a data frame with identifiers for row and column
#'
#' This is meant to be a lightweight replacement for \code{reshape2::melt} that
#' will only work with matrices.
#'
#' @param x A matrix
#'
#' @return A data frame with the number of rows equal to the number of elements
#'   in \code{x}.
#' @export
#'
#' @examples
#' x <- matrix(rnorm(9),3,3)
#' melt_mat(x)
melt_mat <- function(x) {
  out <- data.frame(row = c(row(x)), col = c(col(x)), value = c(x))
  return(out)
}

#' Lightweight tile plot function
#'
#' This function has no dependencies outside of base R, but provides a
#'   reasonable approximation to the functionality of
#'   \code{ggplot2::geom_raster}.
#'
#' @param tile_df A data frame with three columns: \code{row}, \code{col}, and
#'   \code{value} describing locations within a matrix. See
#'   \code{\link{melt_mat}} for an example.
#' @param col A color palette
#' @param ncols The number of colors for a color palette, if \code{col} is not
#'   provided.
#' @param main Plot title (character)
#' @param zlim a vector with the minimum and maximum limits for the plotted
#'   values.
#' @param na.color The color that should be used to represent \code{NA} values
#'   in the tile plot
#'
#' @return A tile plot done in base R graphics
#' @importFrom grDevices heat.colors
#' @importFrom graphics axis layout par rect text
#' @export
#'
#' @examples
#' x <- matrix(rnorm(50*50),50,50)
#' x_df <- melt_mat(x)
#' tile.plot(x_df)
tile.plot <- function(tile_df, col = NULL, ncols = NULL,
                      main = "", zlim = NULL, na.color  = "grey80") {
  .pardefault <- par()
  if("matrix" %in% class(tile_df)) tile_df <- melt_mat(tile_df)
  if(!is.null(col) & !is.null(ncols)) {
    warning("Defining ncols based on col.")
    ncols <- length(col)
  }
  if(is.null(col) & is.null(ncols)) {
    ncols <- 100
    col <- heat.colors(ncols)
  } else {
    if(is.null(ncols)) {
      ncols <- length(col)
    } else {
      col <- heat.colors(ncols)
    }
  }
  if(is.null(zlim)) zlim <- c(min(tile_df$value, na.rm = T), max(tile_df$value, na.rm = T))
  legend_ticks <- round(seq(zlim[1],zlim[2], length.out = 6),
                        digits = 3)
  prob_breaks <- seq(0,1,length.out = ncols)
  pb_diff <- prob_breaks[2] - prob_breaks[1]
  col_quants <- quantile(tile_df$value,na.rm = T,probs = prob_breaks)
  color_breaks <- quantile(seq(zlim[1],zlim[2],length.out = ncols),probs = prob_breaks)
  tile_cols <- vector("numeric",length(tile_df$value))
  tile_cols[tile_df$value >= color_breaks[ncols]] <- col[ncols]
  for(q in rev(seq(ncols))) {
    tile_cols[tile_df$value <= color_breaks[q]] <- col[q]
  }
  tile_cols[tile_cols == "0"] <- na.color
  rows <- max(tile_df$row)
  cols <- max(tile_df$col)
  cb_prime <- min(diff(color_breaks))
  par(mfrow = c(1,2), mar = c(1,1,2,1))
  layout(mat = matrix(c(1,2),nrow = 1, ncol = 2),widths = c(1.7,0.3))
  plot(c(0,rows), c(0,cols), type = 'n', xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", main = main)
  rect(xleft = tile_df$row - 1,ybottom = tile_df$col - 1,
       xright = tile_df$row, ytop = tile_df$col, col = tile_cols,
       border = NA)
  par(mar=c(1,1,2,4))
  plot(c(0,1),zlim, type = "n", xaxt = "n", yaxt = "n", xlab = "",
       ylab = "", bty = "n")
  rect(xleft = 0,
       ybottom = color_breaks,
       xright = 1,
       ytop = color_breaks + cb_prime,
       col = col, border = NA)
  axis(side = 4,at = legend_ticks,
       labels = rep("",6), srt = 45, tck = 0.5)
  text(x = 1, adj = c(-1,0), pos = 4, y = legend_ticks,
       labels = legend_ticks, srt = 0, xpd = NA)
  par(mfrow = .pardefault$mfrow, mar = .pardefault$mar)
}
