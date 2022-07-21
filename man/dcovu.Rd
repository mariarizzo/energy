\name{Unbiased distance covariance}
\alias{bcdcor}
\alias{dcovU}
\title{Unbiased dcov and bias-corrected dcor statistics}
\description{
 These functions compute unbiased estimators of squared distance
 covariance and a bias-corrected estimator of
 (squared) distance correlation.
}
\usage{
bcdcor(x, y)
dcovU(x, y)
}
\arguments{
  \item{x}{ data or dist object of first sample}
  \item{y}{ data or dist object of second sample}
}
\details{
 The unbiased (squared) dcov is inner product definition of
 dCov, in the Hilbert space of U-centered distance matrices.

 The sample sizes (number of rows) of the two samples must
 agree, and samples must not contain missing values. 
 
Argument types supported are 
numeric data matrix, data.frame, or tibble, with observations in rows;
numeric vector; ordered or unordered factors. In case of unordered factors
a 0-1 distance matrix is computed.
}
\value{
\code{dcovU} returns the unbiased estimator of squared dcov.
\code{bcdcor} returns a bias-corrected estimator of squared dcor.
}
\note{
Unbiased distance covariance (SR2014) corresponds to the biased
(original) \eqn{\mathrm{dCov^2}}{dCov^2}. Since \code{dcovU} is an
unbiased statistic, it is signed and we do not take the square root.
For the original distance covariance test of independence (SRB2007,
SR2009), the distance covariance test statistic is the V-statistic
\eqn{\mathrm{n\, dCov^2} = n \mathcal{V}_n^2}{n V_n^2} (not dCov).
Similarly, \code{bcdcor} is bias-corrected, so we do not take the
square root as with dCor.
}
\references{
 Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities.
 \emph{Annals of Statistics}, Vol. 42 No. 6, 2382-2412.

 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \doi{10.1214/009053607000000505}

 Szekely, G.J. and Rizzo, M.L. (2009),
 Brownian Distance Covariance,
 \emph{Annals of Applied Statistics},
 Vol. 3, No. 4, 1236-1265.
 \cr \doi{10.1214/09-AOAS312}
  }
  \author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- iris[1:50, 1:4]
 y <- iris[51:100, 1:4]
 dcovU(x, y)
 bcdcor(x, y)
}
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
