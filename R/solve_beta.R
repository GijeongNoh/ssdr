#' solve beta using lasso estimation
#'
#' @param x The name of predictor matrix
#' @param y The name of response matrix
#' @param G2 The square root of covariance x
#' @param lambda1 shrinkage parameters
#' @param lambda2 shrinkage parameters
#'
#' @import glmnet
#'
#' @export
#'
solve_beta<-function(x, y, G2, lambda1, lambda2)
{
  # transform data to an L1 problem
  x.star<-rbind(x, sqrt(lambda2) * G2)
  y.star<-c(y, rep(0, ncol(x)))

  fit2 <- glmnet::glmnet(x.star, y.star, family="gaussian")
  beta.est <- coef(fit2, s=lambda1)[-1,]

  # return
  return(beta.est)
}
