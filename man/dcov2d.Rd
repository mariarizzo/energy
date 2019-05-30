\name{Fast bivariate dcor and dcov}
\alias{dcor2d}
\alias{dcov2d}
\title{Fast dCor and dCov for bivariate data only}
\description{
For bivariate data only, these are fast O(n log n) implementations of distance
correlation and distance covariance statistics. The U-statistic for dcov^2 is unbiased; 
the V-statistic is the original definition in SRB 2007. These algorithms do not
store the distance matrices, so they are suitable for large samples. 
}
\usage{
dcor2d(x, y, type = c("V", "U"))
dcov2d(x, y, type = c("V", "U"), all.stats = FALSE)
}
\arguments{
  \item{x}{ numeric vector}
  \item{y}{ numeric vector}
  \item{type}{ "V" or "U", for V- or U-statistics}
  \item{all.stats}{ logical}
}
\details{
The unbiased (squared) dcov is documented in \code{dcovU}, for multivariate data in arbitrary, not necessarily equal dimensions. \code{dcov2d} and \code{dcor2d} provide a faster O(n log n) algorithm for bivariate (x, y) only (X and Y are real-valued random vectors). The O(n log n) algorithm was proposed by Huo and Szekely (2016). The algorithm is faster above a certain sample size n. It does not store the distance matrix so the sample size can be very large. 
}
\value{
By default, \code{dcor2d} returns the V-statistic dCor_n^2(x, y), and if type="U", it returns a bias-corrected estimator of squared dcor.

By default, \code{dcov2} returns the V-statistic V_n^2 = dCov_n^2(x, y), and if type="U", it returns the U-statistic, unbiased for V^2(X, Y). The argument all.stats=TRUE is used internally when the function is called from \code{dcor2}. 

For \code{dcov2d} and \code{dcor2d}, direct computation using the C++ function \code{dcovU_stats} may be somewhat faster on small samples, depending on the platform.

\code{dcor2d} and \code{dcov2d} do not store the distance matrices so these functions are helpful when sample size is large, the data is bivariate, and we simply require the statistics. There is not an efficient way to do the nonparametric test by permutations without storing distances. For a test of independence on moderate size samples, use \code{dcov.test} or \code{dcor.test}.
}
\note{
Unbiased distance covariance (SR2014) is equivalent to the U-statistic 
estimator of \eqn{\mathrm{dCov^2}}{dCov^2}. Since \code{dcovU} is an
unbiased statistic, it can be negative and its square root would be 
complex so the square root of the U-statistic is not applied.
For the original distance covariance test of independence (SRB2007,
SR2009), the test statistic was the V-statistic \eqn{\mathrm{n\, dCov_n^2} = n \mathcal{V}_n^2}{n V_n^2}.
Similarly, \code{bcdcor} is bias-corrected, so we do not take the
square root as with \eqn{dCor_n^2}{dCor_n^2}.
}
\references{
Huo, X. and Szekely, G.J. (2016). Fast computing for 
distance covariance. Technometrics, 58(4), 435-447.

 Szekely, G.J. and Rizzo, M.L. (2014),
 Partial Distance Correlation with Methods for Dissimilarities.
 \emph{Annals of Statistics}, Vol. 42 No. 6, 2382-2412.

 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \url{http://dx.doi.org/10.1214/009053607000000505}
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\seealso{
 \code{\link{dcov2d.test}} 
 
 \code{\link{dcov}} \code{\link{dcov.test}}  \code{\link{dcor}} \code{\link{dcor.test}} (multivariate statistics and permutation test) 
}
\examples{
  \donttest{
    ## these are equivalent, but 2d is faster for n > 50
    n <- 100
    x <- rnorm(100)
    y <- rnorm(100)
    all.equal(dcov(x, y)^2, dcov2d(x, y), check.attributes = FALSE)
    all.equal(bcdcor(x, y), dcor2d(x, y, "U"), check.attributes = FALSE)
    }
}
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
