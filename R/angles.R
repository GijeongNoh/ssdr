#' angle between two spaces
#'
#' @param B1 The name of first space
#' @param B2 The name of second space
#'
#' @export
#'
angles<-function(B1, B2)
{
  if(!is.matrix(B1)) B1<-as.matrix(B1)
  if(!is.matrix(B2)) B2<-as.matrix(B2)

  if(ncol(B1) >= ncol(B2)) {
    B<-B1; B.hat<-B2
  } else {
    B<-B2; B.hat<-B1
  }

  P1<-B %*% solve(t(B) %*% B) %*% t(B)
  if(ncol(B.hat) == 1) {
    nume<-as.vector(t(B.hat) %*% P1 %*% B.hat)
    deno<-as.vector(t(B.hat) %*% B.hat)
    ratio<-nume / deno
  } else {
    BtB<-t(B.hat) %*% B.hat
    ei<-eigen(BtB)
    BtB2<-ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
    M<-BtB2 %*% t(B.hat) %*% P1 %*% B.hat %*% BtB2
    ratio<-abs(eigen(M)$values[nrow(M)])
  }
  ans<-acos(sqrt(ratio))/pi * 180
  if(ans > 90) ans<-180 - ans
  return(ans)
}
