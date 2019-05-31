\name{dcov2d}
\alias{dcor2d}
\alias{dcov2d}
\alias{dcov2d.test}
\alias{dcor2d.test}
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
dcov2d.test(x, y, m=250)
dcor2d.test(x, y, m=250)
}
\arguments{
  \item{x}{ numeric vector}
  \item{y}{ numeric vector}
  \item{type}{ "V" or "U", for V- or U-statistics}
  \item{all.stats}{ logical}
  \item{m}{ number of eigenvalues}
}
\details{
The unbiased (squared) dcov is documented in \code{dcovU}, for multivariate data in arbitrary, not necessarily equal dimensions. \code{dcov2d} and \code{dcor2d} provide a faster O(n log n) algorithm for bivariate (x, y) only (X and Y are real-valued random vectors). The O(n log n) algorithm was proposed by Huo and Szekely (2016). The algorithm is faster above a certain sample size n. It does not store the distance matrix so the sample size can be very large. An O(n log n) test based on the asymptotic distribution of \eqn{n V_n}{n V_n} is provided by \code{dcov2d.test} and \code{dcor2d.test}.

For details about inference with distance covariance see \code{dcov}, 
\code{dcor}, \code{dcov.test} and \code{dcor.test}. 

The argument \code{m} specifies the number of eigenvalues estimated for computing the asymptotic null distribution of the test statistic. Use the permutation test \code{dcov.test} (or \code{dcor.test}) if n < 200. Typically m = 250 is large enough for accurate upper tail probabilities.   
}
\value{
By default, \code{dcov2d} returns the V-statistic \eqn{V_n = dCov_n^2(x, y)}{V_n = dCov_n^2(x, y)}, and if type="U", it returns the U-statistic, unbiased for \eqn{dCov^2(X, Y)}{dCov^2(X,Y)}. The argument all.stats=TRUE is used internally when the function is called from \code{dcor2d}. 

By default, \code{dcor2d} returns \eqn{dCor_n^2(x, y)}{dCor_n^2(x, y)}, and if type="U", it returns a bias-corrected estimator of squared dcor equivalent to \code{bcdcor}.

\code{dcov2d.test} and \code{dcor2d.test} return a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{   estimate}{ V_n, R_n (squared dcor)}
   \item{    p.value}{ approximate p-value of the test}
   \item{  data.name}{ description of data}

These functions do not store the distance matrices so they are helpful when sample size is large and the data is bivariate. 
}
\note{
The U-statistic \eqn{U_n}{U_n} can be negative in the lower tail so 
the square root of the U-statistic is not applied.
Similarly, \code{dcor2d(x, y, "U")} is bias-corrected and can be
negative in the lower tail, so we do not take the
square root. The original definitions of dCov and dCor 
(SRB2007, SR2009) were based on V-statistics, which are non-negative,
and defined using the square root of V-statistics.

It has been suggested that instead of taking the square root of the U-statistic, one could take the root of \eqn{|U_n|}{|U_n|} before applying the sign, but that introduces more bias than the original dCor, and should never be used. 
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\seealso{
 \code{\link{dcov}} \code{\link{dcov.test}}  \code{\link{dcor}} \code{\link{dcor.test}} (multivariate statistics and permutation test) 
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
\examples{
  \donttest{
    ## these are equivalent, but 2d is faster for n > 50
    n <- 100
    x <- rnorm(100)
    y <- rnorm(100)
    all.equal(dcov(x, y)^2, dcov2d(x, y), check.attributes = FALSE)
    all.equal(bcdcor(x, y), dcor2d(x, y, "U"), check.attributes = FALSE)

    x <- rlnorm(400)
    y <- rexp(400)
    dcov.test(x, y, R=199)  #permutation test
    dcov2d.test(x, y)       #large sample test
    dcor.test(x, y, R=199)
    dcor2d.test(x, y)
    }  
}
\keyword{ htest }
\keyword{ nonparametric }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
