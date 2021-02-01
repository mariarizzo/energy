\name{Poisson GOF Test}
\alias{poissonGOF.test}
\alias{poisson.mtest}
\alias{poisson.m}
\alias{poisson.e}
\title{ Goodness-of-Fit Test for Poisson Distribution}
\description{
 Performs the mean distance goodness-of-fit test or the energy goodness-of-fit test of Poisson distribution with unknown parameter.
}
\usage{
poissonGOF.test(x, R, method)
poisson.mtest(x, R)
poisson.m(x)
poisson.e(x)
}
\arguments{
  \item{x}{ vector of nonnegative integers, the sample data }
  \item{R}{ number of bootstrap replicates }
  \item{method}{ character: "M" or "E"}
}
\details{
If \code{method} is "M" (default) the Poisson M-test is applied; otherwise the "E" or energy test is applied. Both tests are based on expected distances between observations.

 The mean distance test of Poissonity was proposed and implemented by
 Szekely and Rizzo (2004). The test is based on the result that the sequence
 of expected values E|X-j|, j=0,1,2,... characterizes the distribution of
 the random  variable X. As an application of this characterization one can
 get an estimator \eqn{\hat F(j)} of the CDF. The test statistic
 (see \code{\link{poisson.m}}) is a Cramer-von Mises type of distance, with
 M-estimates replacing the usual EDF estimates of the CDF:
  \deqn{M_n = n\sum_{j=0}^\infty (\hat F(j) - F(j\;; \hat \lambda))^2
 f(j\;; \hat \lambda).}{M_n = n sum [j>=0] (\hat F(j) - F(j; \hat \lambda))^2
 f(j; \hat \lambda).} The test is implemented by parametric bootstrap with
 \code{R} replicates.
 
An energy goodness-of-fit test is based on the test statistic
\deqn{Q_n = n (\frac{2}{n} \sum_{i=1}^n E|x_i - X| - E|X-X'| - \frac{1}{n^2} \sum_{i,j=1}^n |x_i - x_j|,
}{Q_n = n((2/n) sum[1:n] E|x_i-X| - E|X-X'| - (1/n^2) sum[1:n,1:n]
 |x_i-x_j|),}
where X and X' are iid with the hypothesized null distribution. For a test of H: X ~ Poisson(\eqn{\lambda}), we can express E|X-X'| in terms of Bessel functions, and E|x_i - X| in terms of the CDF of Poisson(\eqn{\lambda}).
}

\value{
The functions \code{poisson.m} and \code{poisson.e} return the test statistics. The function
\code{poissonGOF.test} returns a list with class \code{htest} containing
 \item{method}{Description of test}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{description of data}
 \item{estimate}{sample mean}
}
\references{
Szekely, G. J. and Rizzo, M. L. (2004) Mean Distance Test of Poisson Distribution, \emph{Statistics and Probability Letters},
67/3, 241-247. \url{http://dx.doi.org/10.1016/j.spl.2004.01.005}.

Szekely, G. J. and Rizzo, M. L. (2005) A New Test for
 Multivariate Normality, \emph{Journal of Multivariate Analysis},
 93/1, 58-80,
 \url{http://dx.doi.org/10.1016/j.jmva.2003.12.002}.
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- rpois(50, 2)
 poisson.m(x)
 poisson.e(x)
 \donttest{
   set.seed(2021)
   poissonGOF.test(x, R = 199, method="M")
   set.seed(2021)
   poissonGOF.test(x, R = 199, method="E")
 }
}
\keyword{ htest }
\keyword{ energy }