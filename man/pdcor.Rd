\name{pdcor}
\alias{pdcor}
\alias{pdcov}
\alias{pdcor.test}
\alias{pdcov.test}
\title{
Partial distance correlation and covariance
}
\description{Partial distance correlation pdcor, pdcov, and tests.}
\usage{
  pdcov.test(x, y, z, R)
  pdcor.test(x, y, z, R)
  pdcor(x, y, z)
  pdcov(x, y, z)
}
\arguments{
\item{x}{ data or dist object of first sample}
\item{y}{ data or dist object of second sample}
\item{z}{ data or dist object of third sample}
\item{R}{ replicates for permutation test}
}
\details{
\code{pdcor(x, y, z)} and \code{pdcov(x, y, z)} compute the partial distance
correlation and partial distance covariance, respectively,
of x and y removing z.

A test for zero partial distance correlation (or zero partial distance covariance) is implemented in \code{pdcor.test}, and \code{pdcov.test}.

Argument types supported are numeric data matrix, data.frame, tibble, numeric vector, class "dist" object, or factor. For unordered factors a 0-1 distance matrix is computed. 
}
\value{
Each test returns an object of class \code{htest}.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\references{
Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities.
 \emph{Annals of Statistics}, Vol. 42 No. 6, 2382-2412.
}
\examples{
  n = 30
  R <- 199

  ## mutually independent standard normal vectors
  x <- rnorm(n)
  y <- rnorm(n)
  z <- rnorm(n)

  pdcor(x, y, z)
  pdcov(x, y, z)
  set.seed(1)
  pdcov.test(x, y, z, R=R)
  set.seed(1)
  pdcor.test(x, y, z, R=R)

\donttest{
  if (require(MASS)) {
    p = 4
    mu <- rep(0, p)
    Sigma <- diag(p)
  
    ## linear dependence
    y <- mvrnorm(n, mu, Sigma) + x
    print(pdcov.test(x, y, z, R=R))
  
    ## non-linear dependence
    y <- mvrnorm(n, mu, Sigma) * x
    print(pdcov.test(x, y, z, R=R))
    }
  }
}
\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
