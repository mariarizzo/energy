\name{Distance Matrix}
\alias{is.dmatrix}
\alias{calc_dist}
\title{ Distance Matrices }
\description{
Utilities for working with distance matrices. 
\code{is.dmatrix} is a utility that checks whether the argument is a distance or dissimilarity matrix; is it square symmetric, non-negative, with zero diagonal? \code{calc_dist} computes a distance matrix directly from a data matrix. 
}
\usage{
is.dmatrix(x, tol = 100 * .Machine$double.eps)
calc_dist(x)
}
\arguments{
  \item{x}{ numeric matrix}
  \item{tol}{ tolerance for checking required conditions}
} 
\details{
Energy functions work with the distance matrices of samples. The \code{is.dmatrix} function is used internally when converting arguments to distance matrices. The default \code{tol} is the same as default tolerance of \code{isSymmetric}.

\code{calc_dist} is an exported Rcpp function that returns a Euclidean distance matrix from the input data matrix. 
}
\value{
\code{is.dmatrix} returns TRUE if (within tolerance) \code{x} is a distance/dissimilarity matrix; otherwise FALSE. It will return FALSE if \code{x} is a class \code{dist} object.

\code{calc_dist} returns the Euclidean distance matrix for the data matrix \code{x}, which has observations in rows. 
}
\note{
In practice, if \code{dist(x)} is not yet computed, \code{calc_dist(x)} will be faster than \code{as.matrix(dist(x))}. 

On working with non-Euclidean dissimilarities, see the references.
}
\examples{
x <- matrix(rnorm(20), 10, 2)
D <- calc_dist(x)
is.dmatrix(D)
is.dmatrix(cov(x))
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu}
}
\references{
Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities.
 \emph{Annals of Statistics}, Vol. 42 No. 6, 2382-2412.
}
