\name{dcovU_stats}
\alias{dcovU_stats}
\title{Unbiased distance covariance statistics}
\description{
 This function computes unbiased estimators of squared distance
 covariance, distance variance, and a bias-corrected estimator of
 (squared) distance correlation.
}
\usage{
dcovU_stats(Dx, Dy)
}
\arguments{
  \item{Dx}{ distance matrix of first sample}
  \item{Dy}{ distance matrix of second sample}
}
\details{
 The unbiased (squared) dcov is inner product definition of
 dCov, in the Hilbert space of U-centered distance matrices.

 The sample sizes (number of rows) of the two samples must
 agree, and samples must not contain missing values. The
 arguments must be square symmetric matrices. 
}
\value{
\code{dcovU_stats} returns a vector of the components of bias-corrected
dcor: [dCovU, bcdcor, dVarXU, dVarYU].
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
 Dx <- as.matrix(dist(x))
 Dy <- as.matrix(dist(y))
 dcovU_stats(Dx, Dy)
 }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
