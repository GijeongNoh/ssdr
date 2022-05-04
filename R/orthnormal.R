#' Gram-Schmidt orthonormalization
#'
#' @param X to gram-schmidt orthonormalization
#'
#' @export
orthnormal<-function(X)
{
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)

  W<-NULL
  if(p > 1) {
    W<-cbind(W, X[,1])
    for(k in 2:p) {
      gw<-rep(0, n)
      for(i in 1:(k-1)) {
        gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
        gw<-gw + gki * W[,i]
      }
      W<-cbind(W, X[,k] - gw)
    }
  } else {
    W<-cbind(W, X[,1])
  }

  W<-apply(W, 2, norm)
  W
}
