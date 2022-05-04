#' distance between subspaces
#'
#' @param v1 The name of first subspace
#' @param v2 The name of second subspace
#'
#' @export
dist=function(v1,v2){
  v1=as.matrix(v1);v2=as.matrix(v2)
  if(dim(v1)[1]>1){
    p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
    p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
  }
  if(dim(v1)[1]==1){
    p1=v1%*%t(v1)/c(t(v1)%*%v1)
    p2=v2%*%t(v2)/c(t(v2)%*%v2)}
  d <- sqrt(sum((p1-p2)*(p1-p2)))
  return(d)
}
