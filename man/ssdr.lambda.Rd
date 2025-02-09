% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssdr.lambda.R
\name{ssdr.lambda}
\alias{ssdr.lambda}
\title{Sparse sufficient dimension reduction with shrinkage parameter lambda}
\usage{
ssdr.lambda(
  X,
  y,
  method = c("pc", "sir", "save", "phdres", "dr"),
  d,
  nslices,
  lambda1,
  lambda2,
  max.iter,
  eps.conv = 0.001
)
}
\arguments{
\item{X}{The name of predictor matrix}

\item{y}{The name of response matrix}

\item{method}{sufficient dimension reduction method}

\item{d}{The number of dimensions to reduce}

\item{nslices}{The number of slices}

\item{lambda1}{shrinkage parameters}

\item{lambda2}{shrinkage parameters}

\item{max.iter}{The number of max iterations}

\item{eps.conv}{Epsilon convergence}
}
\value{
The output from sparse sufficient dimension reduction
}
\description{
Sparse sufficient dimension reduction with shrinkage parameter lambda
}
