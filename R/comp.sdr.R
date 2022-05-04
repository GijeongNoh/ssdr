#' Compute with sufficient dimension reduction method
#'
#' @param X The name of predictor matrix
#' @param y The name of response matrix
#' @param method sufficient dimension reduction method
#' @param d The number of dimensions to reduce
#' @param nslices The number of slices
#'
#' @import stats
#' @import dr
#'
#' @return The output from sufficient dimension reduction
#' @export
comp.sdr<-function(X, y, method, d, nslices)
{
  out.m<-switch(method,
                pc     = comp.sdr.pc(X, y, d),
                sir    = comp.sdr.sir(X, y, d, nslices),
                save   = comp.sdr.save(X, y, d, nslices),
                phdres = comp.sdr.phdres(X, y, d),
                dr     = comp.sdr.dr(X, y, d, nslices)
  )

  ans<-list(m=out.m$m, G=out.m$G, beta.sdr=out.m$beta.sdr)
  return(ans)
}



comp.sdr.pc<-function(X, y, d)
{
  # m matrix
  Xc<-apply(X, 2, center_vec)
  m<-mat.sqrt(t(Xc) %*% Xc)

  # G matrix
  G<-diag(1, ncol(X))

  # pc
  v<-eigen(cov.x(X))$vectors[, 1:d]
  if(d == 1) v<-matrix(v, ncol=1)

  # return
  ans<-list(m=m, G=G, beta.sdr=v)
  return(ans)
}



comp.sdr.sir<-function(X, y, d, nslices)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2

  # slice y
  sy<-dr.slices(y, nslices)
  nslices<-sy$nslices

  # compute sdr kernel matrix
  M.sir.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:nslices) {
    Z.s<-Z[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
    Z.sm<-as.vector(apply(Z.s, 2, mean))
    M.sir.z<-M.sir.z + (sy$slice.sizes[s]/n) * Z.sm %*% t(Z.sm) #sirmat?? ?????? candidate matrix
  }
  M.sir<-Sigma.x2 %*% M.sir.z %*% Sigma.x2
  #m<-mat.sqrt(n*M.sir)
  m<-mat.sqrt(M.sir) #M^(1/2) -> mi?? ???ϱ? ��??(ssdr�� ��?ؼ?)

  # compute sdr estimate w/o shrinkage(?׳? sdr�� ��?ؼ?)
  v<-eigen(M.sir.z)$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)

  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}



comp.sdr.save<-function(X, y, d, nslices)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2

  # slice y
  sy<-dr.slices(y, nslices)
  nslices<-sy$nslices

  # compute sdr kernel matrix
  M.save.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:nslices) {
    Z.s<-Z[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
    iVz<-diag(1, p) - cov.x(Z.s)
    M.save.z<-M.save.z + (sy$slice.sizes[s]/n) * iVz %*% iVz
  }
  M.save<-Sigma.x2 %*% M.save.z %*% Sigma.x2
  m<-mat.sqrt(M.save)

  # compute sdr estimate w/o shrinkage
  v<-eigen(M.save.z)$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)

  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}



comp.sdr.phdres<-function(X, y, d)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov.x(X)
  Sigma.x2<-mat.sqrt(Sigma.x)
  Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

  # standardize X
  Z<-apply(X, 2, center_vec) %*% Sigma.x.inv2

  # residual
  e<-resid(lm(as.vector(y)~Z))

  # compute sdr kernel matrix
  M.phd.z<-matrix(0, nrow=p, ncol=p)
  for(s in 1:n) {
    M.phd.z<-M.phd.z + e[s] * Z[s,] %*% t(Z[s,])
  }
  M.phd.z<-M.phd.z / n
  M.phd<-Sigma.x2 %*% (M.phd.z %*% t(M.phd.z)) %*% Sigma.x2
  m<-mat.sqrt(M.phd)

  # compute sdr estimate w/o shrinkage
  v<-eigen(M.phd.z %*% t(M.phd.z))$vectors[,1:d]
  if(d == 1) v<-matrix(v, ncol=1)
  beta.sdr<-Sigma.x.inv2 %*% v
  beta.sdr<-apply(beta.sdr, 2, norm)

  # return
  ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
  return(ans)
}


comp.sdr.dr <- function(X,y,d,nslices){
  p=ncol(X);n=nrow(X)
  Sigma.x=cov.x(X)
  Sigma.x2=mat.sqrt(Sigma.x)
  signrt=mat.sqrt.inv(Sigma.x)
  xc=t(t(X)-apply(X,2,mean))
  xst=xc%*%signrt
  ydis=discretize(y,nslices)
  ylabel=unique(ydis)
  prob=numeric()
  for(i in 1:nslices) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,nslices));exy=numeric()
  for(i in 1:nslices) {
    vxy[,,i]=var(xst[ydis==ylabel[i],])
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:nslices){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  M.dr = Sigma.x2 %*% out %*% Sigma.x2
  m = mat.sqrt(M.dr)
  ans <- list(m=m, G=Sigma.x, beta.sdr =signrt%*%eigen(out)$vectors[,1:d])
  return(ans)
}
