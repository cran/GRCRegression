\name{link.log}
\alias{link.log}
\title{The log link}
\description{The log link function}
\details{
  This implements the link function \eqn{g(x)=\log(x)}{g(x) = log(x)}
  on the domain \eqn{(0, \infty)}{(0, Inf)}.
}
\usage{
link.log()
}
\arguments{No arguments}
\value{
  \item{class}{A character \code{"0Inf"} denotes the domain \eqn{(0,\infty)}{(0, Inf)}.}
  \item{g}{The link function \eqn{g}{g}.}
  \item{gInv}{The inverse of the link function.}
  \item{D.gInv}{The derivative of the inverse link function.}
  \item{DD.gInv}{The second-order derivative of the link function.}
}

\examples{
a <- link.log()
}
