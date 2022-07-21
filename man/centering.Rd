\name{centering distance matrices}
\alias{Ucenter}
\alias{Dcenter}
\alias{U_center}
\alias{D_center}
\title{ Double centering and U-centering }
\description{
 Stand-alone double centering and U-centering functions
 that are applied in unbiased distance covariance, bias
 corrected distance correlation, and partial distance correlation.
 }
\usage{
Dcenter(x)
Ucenter(x)
U_center(Dx)
D_center(Dx)
}
\arguments{
  \item{x}{ dist object or data matrix}
  \item{Dx}{ distance or dissimilarity matrix}
}
\details{
 In \code{Dcenter} and \code{Ucenter}, \code{x} must be
 a \code{dist} object or a data matrix. Both functions return
 a doubly centered distance matrix.

 Note that \code{pdcor}, etc. functions include the
 centering operations (in C), so that these stand alone versions
 of centering functions are not needed except in case one
 wants to compute just a double-centered or U-centered matrix.

 \code{U_center} is the Rcpp export of the cpp function.
 \code{D_center} is the Rcpp export of the cpp function.
}
\value{
All functions return a square symmetric matrix.

 \code{Dcenter} returns a matrix
 \deqn{A_{ij}=a_{ij} - \bar a_{i.} - \bar a_{.j} + \bar a_{..}}
 as in classical multidimensional scaling. \code{Ucenter}
 returns a matrix
 \deqn{\tilde A_{ij}=a_{ij} - \frac{a_{i.}}{n-2}
 - \frac{a_{.j}}{n-2} + \frac{a_{..}}{(n-1)(n-2)},\quad i \neq j,}
 with zero diagonal,
 and this is the double centering applied in \code{pdcov} and
 \code{pdcor} as well as the unbiased dCov and bias corrected
 dCor statistics.
}
\note{
The c++ versions \code{D_center} and \code{U_center} should typically
be faster. R versions are retained for historical reasons.
}
\references{
  Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities,
 \emph{Annals of Statistics}, Vol. 42, No. 6, pp. 2382-2412.
 \cr \url{https://projecteuclid.org/euclid.aos/1413810731}
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- iris[1:10, 1:4]
 dx <- dist(x)
 Dx <- as.matrix(dx)
 M <- U_center(Dx)

 all.equal(M, U_center(M))     #idempotence
 all.equal(M, D_center(M))     #invariance
}
\keyword{ multivariate }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
