#' Power of matrix
#'
#' @param a The name of matrix
#' @param alpha The order you want
#'
#' @export
#'
matpower=function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}
