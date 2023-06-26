\name{summary}
\alias{summary.GRCglm.P}
\title{Summary of estimates.}
\description{A summary of estimates for the class \code{GRCglm.P}.}
\usage{
\method{summary}{GRCglm.P}(object, level = 0.95, \dots)
}
\arguments{
  \item{object}{An object of class ``\code{GRCglm.P}''.}
  \item{level}{Level of confidence.}
  \item{\dots}{Reserved parameter for matching other summary methods.}
}
\value{
  \item{FisherInfo}{The observed Fisher information matrix.}
  \item{stdErr}{Estimated standard error.}
  \item{ci}{Confidence intervals of the estimated coefficients.}
  \item{WOGRCStdErr}{Estimated standard error for comparision, with Fisher 
      information matrix computed using exact enumeration of counts.}
  \item{zValue}{Z-score of the estimated coefficients.}
  \item{sigLevel}{P-value of the estimated coefficients.}
}

\examples{
set.seed(123)
tp <- genData.P(beta = c(0.5, -1, 1), data.size = 120, scheme = c(0:3, 5, 8),
  scope.lambda = c(1, 10))
a <- GRCglm.P(y = tp$y, x1 = tp$x, scheme = c(0:3, 5, 8))
b <- summary(a)
}
