#' choose the order of beta
#'
#' @param beta the result beta of sparse sufficient dimension reduction
#'
#' @export
beta.order <- function(beta){
  new.beta <- beta
  if (abs(beta[1,1]) < abs(beta[1,2])){
    new.beta <- matrix(c(beta[,2], beta[,1]), ncol=2)
  }
  return(new.beta)
}
