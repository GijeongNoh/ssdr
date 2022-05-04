#' covariance matrix
#'
#' @param X The name of matrix to find covariance matrix
#'
#' @export
cov.x<-function(X)
{
  Xc<-apply(X, 2, center_vec)
  t(Xc) %*% Xc / nrow(Xc)
}
