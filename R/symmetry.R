#' symmtrize a matrix
#'
#' @param a The name of matrix for symmtrize
#'
#' @export
symmetry = function(a){
  return((a + t(a))/2)}
