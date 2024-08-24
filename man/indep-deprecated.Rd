\name{indep.test}
\alias{indep.test}
\title{ Energy-tests of Independence}
\description{
 Computes a multivariate nonparametric test of independence.
 The default method implements the distance covariance test
 \code{\link{dcov.test}}.
 }
\usage{
indep.test(x, y, method = c("dcov","mvI"), index = 1, R)
}
\arguments{
  \item{x}{ matrix: first sample, observations in rows}
  \item{y}{ matrix: second sample, observations in rows}
  \item{method}{ a character string giving the name of the test}
  \item{index}{ exponent on Euclidean distances}
  \item{R}{ number of replicates}
}
\details{
 \code{indep.test} with the default \code{method = "dcov"} computes
 the distance
 covariance test of independence. \code{index} is an exponent on
 the Euclidean distances. Valid choices for \code{index} are in (0,2],
 with default value 1 (Euclidean distance). The arguments are passed
 to the \code{dcov.test} function. See the help topic \code{\link{dcov.test}} for
 the description and documentation and also see the references below.

 \code{indep.test} with \code{method = "mvI"}
 computes the coefficient \eqn{\mathcal I_n}{I_n} and performs a nonparametric
 \eqn{\mathcal E}{E}-test of independence. The arguments are passed to
 \code{mvI.test}. The
 \code{index} argument is ignored (\code{index = 1} is applied).
 See the help topic \code{\link{mvI.test}} and also
 see the reference (2006) below for details.

 The test decision is obtained via
 bootstrap, with \code{R} replicates.
 The sample sizes (number of rows) of the two samples must agree, and
 samples must not contain missing values.

 These energy tests of independence are based on related theoretical
 results, but different test statistics.
 The \code{dcov} method is faster than \code{mvI} method by
 approximately a factor of O(n).
}
\value{
 \code{indep.test} returns a list with class
 \code{htest} containing
 \item{     method}{description of test}
 \item{  statistic}{observed value of the
 test statistic \eqn{n \mathcal V_n^2}{n V_n^2}
 or \eqn{n \mathcal I_n^2}{n I_n^2}}
    \item{  estimate}{ \eqn{\mathcal V_n}{V_n} or \eqn{\mathcal I_n}{I_n}}
    \item{  estimates}{ a vector [dCov(x,y), dCor(x,y), dVar(x), dVar(y)]
    (method dcov)}
   \item{ replicates}{ replicates of the test statistic}
   \item{    p.value}{approximate p-value of the test}
 \item{  data.name}{description of data}
 }
\note{As of energy-1.1-0,
\code{indep.etest} is deprecated and replaced by \code{indep.test}, which
has methods for two different energy tests of independence.  \code{indep.test} applies
the distance covariance test (see \code{dcov.test}) by default (\code{method = "dcov"}).
The original \code{indep.etest} applied the independence coefficient
\eqn{\mathcal I_n}{I_n}, which is now obtained by \code{method = "mvI"}.
 }
 \seealso{
 \code{ \link{dcov.test} }
 \code{ \link{mvI.test} }
 \code{ \link{dcov} }
 \code{ \link{mvI} }
 }
\references{
  Szekely, G.J. and Rizzo, M.L. (2009),
  Brownian Distance Covariance,
 \emph{Annals of Applied Statistics}, Vol. 3 No. 4, pp.
 1236-1265. (Also see discussion and rejoinder.)
 \cr \doi{10.1214/09-AOAS312}

  Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
  Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \doi{10.1214/009053607000000505}

 Bakirov, N.K., Rizzo, M.L., and Szekely, G.J. (2006), A Multivariate
 Nonparametric Test of Independence, \emph{Journal of Multivariate Analysis}
 93/1, 58-80, \cr
 \doi{10.1016/j.jmva.2005.10.005}
 }
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
\donttest{
 ## independent multivariate data
 x <- matrix(rnorm(60), nrow=20, ncol=3)
 y <- matrix(rnorm(40), nrow=20, ncol=2)
 indep.test(x, y, method = "dcov", R = 99)
 indep.test(x, y, method = "mvI", R = 99)

 ## dependent multivariate data
 if (require(MASS)) {
   Sigma <- matrix(c(1, .1, 0, 0 , 1, 0, 0 ,.1, 1), 3, 3)
   x <- mvrnorm(30, c(0, 0, 0), diag(3))
   y <- mvrnorm(30, c(0, 0, 0), Sigma) * x
   indep.test(x, y, R = 99)    #dcov method
   indep.test(x, y, method = "mvI", R = 99)
    }
  }
}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ energy statistics }

