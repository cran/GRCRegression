\name{genData.P}
\alias{genData.P}
\title{Artificial Data Generation}
\description{Generate artificial data for testing and verification.}
\usage{
genData.P(beta, data.size, scheme,
  scope.lambda = c(0.01, 100), link.lambda = link.log, seed = list(no = 1),
    y.lowRatio = -1)
}
\arguments{
  \item{beta}{Coefficients of modified Poisson regression.}
  \item{data.size}{The sample size of data to generate.}
  \item{scheme}{The grouped and right-censored grouping scheme. A sorted vector
    containing starting integers of all groups.}
  \item{scope.lambda}{A vector of length 2 specifying the range of
    \eqn{\lambda}{lambda}.}
  \item{link.lambda}{The link function for \eqn{\lambda}{lambda}.}
  \item{seed}{Random seed (optional).}
  \item{y.lowRatio}{The minimum proportion of groups in the sample
    distribution. This optional proportion is included to deal with imbalanced
    data.}
}
\value{
  \item{x}{The design matrix.}
  \item{y}{The GRC outcome.}
  \item{unused}{An integer roughly describing the number of discarded sample
  observations in the generating process. If this number is too large, the
  generation of artificial data is not efficient and takes too much time.}
}
\examples{
set.seed(123)
tp <- genData.P(beta = c(0.5, -1, 1), data.size = 12, scheme = c(0:3, 5, 8))
}
