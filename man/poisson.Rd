\name{Poisson Tests}
\alias{poisson.tests}
\alias{poisson.e}
\alias{poisson.etest}
\alias{poisson.m}
\alias{poisson.mtest}
\title{ Goodness-of-Fit Tests for Poisson Distribution}
\description{
Performs the mean distance goodness-of-fit test and the energy goodness-of-fit test of Poisson distribution with unknown parameter.
}
\usage{
poisson.e(x)
poisson.m(x)
poisson.etest(x, R)
poisson.mtest(x, R)
poisson.tests(x, R, test="all")
}
\arguments{
  \item{x}{ vector of nonnegative integers, the sample data }
  \item{R}{ number of bootstrap replicates }
  \item{test}{ name of test(s) }
}
\details{
Two distance-based tests of Poissonity are applied in \code{poisson.tests}, "M" and "E". The default is to 
do all tests and return results in a data frame. 
Valid choices for \code{test} are "M", "E", or "all" with 
default "all". 

If "all" tests, all tests are performed by a single parametric bootstrap computing all test statistics on each sample. 

The "M" choice is two tests, one based on a Cramer-von Mises  distance and the other an Anderson-Darling distance. The "E" choice is the energy goodness-of-fit test. 

\code{R} must be a positive integer for a test. If \code{R} is missing or 0, a warning is printed but test statistics are computed (without testing). 

The mean distance test of Poissonity (M-test) is based on the result that the sequence
 of expected values E|X-j|, j=0,1,2,... characterizes the distribution of
 the random  variable X. As an application of this characterization one can
 get an estimator \eqn{\hat F(j)} of the CDF. The test statistic
 (see \code{\link{poisson.m}}) is a Cramer-von Mises type of distance, with
 M-estimates replacing the usual EDF estimates of the CDF:
  \deqn{M_n = n\sum_{j=0}^\infty (\hat F(j) - F(j\;; \hat \lambda))^2
 f(j\;; \hat \lambda).}{M_n = n sum [j>=0] (\hat F(j) - F(j; \hat \lambda))^2
 f(j; \hat \lambda).} 
 
In \code{poisson.tests}, an Anderson-Darling type of weight is also applied when \code{test="M"} or \code{test="all"}.
 
The tests are implemented by parametric bootstrap with
 \code{R} replicates.
 
An energy goodness-of-fit test (E) is based on the test statistic
\deqn{Q_n = n (\frac{2}{n} \sum_{i=1}^n E|x_i - X| - E|X-X'| - \frac{1}{n^2} \sum_{i,j=1}^n |x_i - x_j|,
}{Q_n = n((2/n) sum[1:n] E|x_i-X| - E|X-X'| - (1/n^2) sum[1:n,1:n]
 |x_i-x_j|),}
where X and X' are iid with the hypothesized null distribution. For a test of H: X ~ Poisson(\eqn{\lambda}), we can express E|X-X'| in terms of Bessel functions, and E|x_i - X| in terms of the CDF of Poisson(\eqn{\lambda}).

If test=="all" or not specified, all tests are run with a single parametric bootstrap. \code{poisson.mtest} implements only the Poisson M-test with Cramer-von Mises type distance. \code{poisson.etest} implements only the Poisson energy test.
}
\value{
The functions \code{poisson.m} and \code{poisson.e} return the test statistics. The function
\code{poisson.mtest} or \code{poisson.etest} return an \code{htest} object containing
 \item{method}{Description of test}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{replicates R}
 \item{estimate}{sample mean}
 
\code{poisson.tests} returns "M-CvM test", "M-AD test" and "Energy test" results in a data frame with columns 
 \item{estimate}{sample mean}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{method}{Description of test}
which can be coerced to a \code{tibble}.
}
\note{The running time of the M test is much faster than the E-test.}
\references{
Szekely, G. J. and Rizzo, M. L. (2004) Mean Distance Test of Poisson Distribution, \emph{Statistics and Probability Letters},
67/3, 241-247. \doi{10.1016/j.spl.2004.01.005}.

Szekely, G. J. and Rizzo, M. L. (2005) A New Test for
 Multivariate Normality, \emph{Journal of Multivariate Analysis},
 93/1, 58-80,
 \doi{10.1016/j.jmva.2003.12.002}.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- rpois(50, 2)
 poisson.m(x)
 poisson.e(x)
 \donttest{
 poisson.etest(x, R=199)
 poisson.mtest(x, R=199)
 poisson.tests(x, R=199)
 }
}
\keyword{ htest }
\keyword{ energy }
