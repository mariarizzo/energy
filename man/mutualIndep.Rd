\name{mutual independence}
\alias{mutualIndep.test}
\title{ Energy Test of Mutual Independence}
\description{
 The test statistic is the sum of d-1 bias-corrected squared dcor statistics where the number of variables is d. Implementation is by permuation test. 
}
\usage{
mutualIndep.test(x, R)
}
\arguments{
  \item{x}{ data matrix or data frame}
  \item{R}{ number of permutation replicates}
}
\details{
A population coefficient for mutual independence of d random variables, \eqn{d \geq 2},  is
\deqn{
  \sum_{k=1}^{d-1} \mathcal R^2(X_k, [X_{k+1},\dots,X_d]).
}
which is non-negative and equals zero iff mutual independence holds. 
For example, if d=4 the population coefficient is 
\deqn{
\mathcal R^2(X_1, [X_2,X_3,X_4]) +
\mathcal R^2(X_2, [X_3,X_4]) +
\mathcal R^2(X_3, X_4),
}
A permutation test is implemented based on the corresponding sample coefficient. 
To test mutual independence of \deqn{X_1,\dots,X_d} the test statistic is the sum of the d-1 
statistics (bias-corrected \eqn{dcor^2} statistics):
\deqn{\sum_{k=1}^{d-1} \mathcal R_n^*(X_k, [X_{k+1},\dots,X_d])}. 
}
\value{
\code{mutualIndep.test} returns an object of class \code{power.htest}.
}
\note{
See Szekely and Rizzo (2014) for details on unbiased \eqn{dCov^2} and bias-corrected \eqn{dCor^2} (\code{bcdcor}) statistics.
}
\seealso{
\code{\link{bcdcor}}, \code{\link{dcovU_stats}}
}
\references{
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \doi{10.1214/009053607000000505}
 
Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities.
 \emph{Annals of Statistics}, Vol. 42 No. 6, 2382-2412.
}

\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
x <- matrix(rnorm(100), nrow=20, ncol=5)
mutualIndep.test(x, 199)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

