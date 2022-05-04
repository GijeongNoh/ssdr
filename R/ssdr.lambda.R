#' Sparse sufficient dimension reduction with shrinkage parameter lambda
#'
#' @param X The name of predictor matrix
#' @param y The name of response matrix
#' @param method sufficient dimension reduction method
#' @param d The number of dimensions to reduce
#' @param nslices The number of slices
#' @param lambda1 shrinkage parameters
#' @param lambda2 shrinkage parameters
#' @param max.iter The number of max iterations
#' @param eps.conv Epsilon convergence
#'
#' @import glmnet
#'
#' @return The output from sparse sufficient dimension reduction
#' @export
ssdr.lambda<-function(X, y, method=c("pc", "sir", "save", "phdres", "dr"), d, nslices, lambda1, lambda2, max.iter, eps.conv=1e-3)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  method<-match.arg(method)

  # compute sdr components
  out.m<-comp.sdr(X, y, method, d, nslices)
  m<-out.m$m
  M<-t(m) %*% m
  G<-out.m$G
  G2<-mat.sqrt(G)
  G2.inv<-mat.sqrt.inv(G2)

  # initial estimate of alpha and beta
  alpha<-out.m$beta.sdr
  beta<-alpha
  for(i in 1:d) {
    ym<-m %*% alpha[,i]
    beta[,i]<-solve_beta(m, ym, G2, lambda1[i], lambda2)
    #fit1 <- glmnet(m, ym, family="gaussian")
    #beta[,i] <- coef(fit1, s=lambda1[i])[-1,]
  }

  # iteration
  iter<-0
  beta.n<-apply(beta, 2, norm)
  diff.conv<-1
  while((iter < max.iter) & (diff.conv[iter+1] > eps.conv)){
    z<-svd(G2.inv %*% M %*% beta)
    alpha<-G2.inv %*% (z$u) %*% t(z$v)
    for(i in 1:d) {
      ym<-m %*% alpha[,i]
      beta[,i]<-solve_beta(m, ym, G, lambda1[i], lambda2)
      #fit1 <- glmnet(m, ym, family="gaussian")
      #beta[,i] <- coef(fit1, s=lambda1[i])[-1,]
    }

    beta.n.new<-apply(beta, 2, norm)
    diff.conv<-c(diff.conv, max(abs(beta.n.new - beta.n)))
    beta.n<-beta.n.new

    iter<-iter + 1
  }

  # compute objective value
  comp1<-m %*% solve(G) - m %*% beta %*% t(beta)
  rss1<-sum(diag(comp1 %*% G %*% t(comp1)))

  z<-svd(G2.inv %*% M %*% beta)
  alpha<-G2.inv %*% (z$u) %*% t(z$v)
  comp2<-m %*% solve(G) - m %*% beta %*% t(alpha)
  rss2<-sum(diag(comp2 %*% G %*% t(comp2)))

  p.e<-sum(as.vector(beta) != 0)

  # normalize beta
  beta<-apply(beta, 2, norm)

  # return
  ans<-list(beta=beta, beta0=out.m$beta.sdr, alpha=alpha, diff.conv=diff.conv, iter=iter, p.e=p.e, rss1=rss1, rss2=rss2)
  return(ans)
}
