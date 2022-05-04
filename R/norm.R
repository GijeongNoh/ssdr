#' normalize a vector
#'
#' @param v The name of the vector for normalization
#'
#' @export
#'
norm<-function(v)
{
  sumv2<-sum(v^2)
  if(sumv2 == 0) sumv2<-1
  v/sqrt(sumv2)
}
