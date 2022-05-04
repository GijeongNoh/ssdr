#' Function to get optimal lambda
#'
#' @param X The name of predictor matrix
#' @param y The name of response matrix
#' @param method sufficient dimension reduction method
#' @param d The number of dimensions to reduce
#' @param nslices The number of slices
#' @param s1.range shrinkage parameters
#' @param s2.range shrinkage parameters
#' @param max.iter The number of max iterations
#' @param eps.conv Epsilon convergence
#'
#' @return Provide optimal lambda
#' @export
ssdr.wrap<-function(X, y, method=c("pc", "sir", "save", "phdres","dr"), d, nslices, s1.range, s2.range, max.iter, eps.conv)
{
  # parameters
  n<-nrow(X)

  # compute criteria for all s1 and s2
  crit.all<-list()
  for(j in 1:length(s2.range)) {
    s2<-s2.range[j]

    crit<-NULL
    for(k in 1:length(s1.range)) {
      s1<-s1.range[k]

      out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)

      aic1<-n*out1$rss1 + 2 * out1$p.e
      bic1<-n*out1$rss1 + log(n) * out1$p.e
      out.crit<-c(aic1, bic1)
      crit<-cbind(crit, out.crit)
    }

    crit.all[[j]]<-crit
  }

  # locate optimal s
  s1.min<-NULL
  ct.min<-NULL
  for(j in 1:length(s2.range)) {
    s1.min.j<-ct.min.j<-NULL
    for(l in 1:length(out.crit)) {
      s1.min.j<-c(s1.min.j, s1.range[order(crit.all[[j]][l,])[1]])
      ct.min.j<-c(ct.min.j, min(crit.all[[j]][l,], na.rm=T))
    }
    s1.min<-rbind(s1.min, s1.min.j)
    ct.min<-rbind(ct.min, ct.min.j)
  }
  rownames(s1.min)<-as.character(s2.range); colnames(s1.min)<-c("aic1", "bic1")
  rownames(ct.min)<-as.character(s2.range); colnames(ct.min)<-c("aic1", "bic1")

  # beta estimate with given optimal s
  beta.est.all<-NULL
  s12.est.all<-NULL
  for(l in 1:length(out.crit)) {
    pos<-order(ct.min[, l])[1]
    s1<-s1.min[pos, l]
    s2<-s2.range[pos]
    out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
    beta.est.all<-cbind(beta.est.all, out1$beta)
    s12.est.all <-cbind(s12.est.all,  c(s1, s2))
  }
  rownames(s12.est.all)<-c("s1", "s2"); colnames(s12.est.all)<-c("aic1", "bic1")

  # return
  ans<-list(beta.est.all=beta.est.all, beta.est0=out1$beta0, s12.est.all=s12.est.all, crit.all=crit.all, s1.min=s1.min, ct.min=ct.min)
  return(ans)
}
