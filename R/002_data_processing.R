#' Format tensor regression data
#'
#' @param y a numeric response of length N
#' @param X an array in which the last margin is of size N
#' @param eta a matrix with N rows
#'
#' @return a list of class 'TR_data' in which dimensions have been verified as
#'   compatible with the tensor regression functions
#' @export
as.TR_data <- function(y, X, eta = NULL) {
  if(length(y) != tail(dim(X),1))
    stop("The last dimension of the array 'X' should be the same size as length(y).")
  if(is.null(eta)) {
    out <- list(y=y,X=X)
  }
  if(!is.null(eta)) {
    eta <- as.matrix(eta)
    if(nrow(eta) != tail(dim(X),1))
      stop("The last dimension of the array 'X' should be the same size as nrow(as.matrix(eta)).")
    out <- list(y=y,X=X,eta=eta)
  }
  class(out) <- "TR_data"
  return(out)
}
