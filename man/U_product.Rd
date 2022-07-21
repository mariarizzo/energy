\name{U_product}
\alias{U_product}
\title{ Inner product in the Hilbert space of U-centered
distance matrices}
\description{
 Stand-alone function to compute the inner product in the
 Hilbert space of U-centered distance matrices, as in the definition of
 partial distance covariance.
}
\usage{
U_product(U, V)
}
\arguments{
  \item{U}{ U-centered distance matrix}
  \item{V}{ U-centered distance matrix}
}
\details{
 Note that \code{pdcor}, etc. functions include the centering and
 projection operations, so that these stand alone versions are not
 needed except in case one wants to check the internal computations.

 Exported from U_product.cpp.
}
\value{
\code{U_product} returns the inner product, a scalar.
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
 y <- iris[11:20, 1:4]
 M1 <- as.matrix(dist(x))
 M2 <- as.matrix(dist(y))
 U <- U_center(M1)
 V <- U_center(M2)

 U_product(U, V)
 dcovU_stats(M1, M2)
 }
\keyword{ multivariate }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
