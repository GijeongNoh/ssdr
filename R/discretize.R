#' discretize
#'
#' @param y The name of y vector
#' @param h The number of slices
#'
#' @export
#'
discretize=function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}
