\name{GRCglm}
\alias{GRCglm.P}
\title{Maximum likelihood estimation for modified Poisson regression of GRC data}
\description{This function implements maximum likelihood estimation of modified
 Poisson regression of grouped and right-censored counts.}
\usage{
  GRCglm.P(y, x1, scheme, link.lambda = link.log, weights = rep(1, nrow(x1)),
    num.intercept = 1, xtol_rel = 1e-08, maxit = 100)
}
\arguments{
  \item{y}{A vector of the GRC outcome.}
  \item{x1}{The design matrix.}
  \item{scheme}{A vector (sorted) of the starting integers of all GRC groups.}
  \item{link.lambda}{The link function for \eqn{\lambda}{lambda}}
  \item{weights}{The weight vector used to consider sampling weights.}
  \item{num.intercept}{Presence of the regression intercept. Can be set to 0 or 1.
    This is also used to calculate the null deviance.}
  \item{xtol_rel, maxit}{The tolerancethreshold and maximum number of iteration of
    the algorithm.}
}
\value{
  \item{coefficients}{The inferred coefficients.}
  \item{beta}{An alias of \code{coefficients}.}
  \item{fitting}{The call.}
  \item{log.likelihood}{Value of the log likelihood function.}
  \item{df.null}{The residual degrees of freedom of the null model.}
  \item{df.residual}{The residual degrees of freedom.}
  \item{null.deviance}{Null deviance.}
  \item{deviance}{The deviance.}
  \item{aic, bic}{The AIC and BIC scores.}
  \item{McFaddenR2, McFaddenAdjR2}{The (rep. Adjusted) McFadden R-square.}
}

\examples{
set.seed(123)
tp <- genData.P(beta = c(0.5, -1, 1), data.size = 120, scheme = c(0:3, 5, 8),
  scope.lambda = c(1, 10))
a <- GRCglm.P(y = tp$y, x1 = tp$x, scheme = c(0:3, 5, 8))
}
