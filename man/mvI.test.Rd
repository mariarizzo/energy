\name{mvI.test}
\alias{mvI.test}
\alias{mvI}
\title{ Energy Statistic Test of Independence}
\description{
 Computes the multivariate nonparametric E-statistic and test of independence
 based on independence coefficient \eqn{\mathcal I_n}{I_n}.}
\usage{
    mvI.test(x, y, R)
    mvI(x, y)
}
\arguments{
  \item{x}{ matrix: first sample, observations in rows}
  \item{y}{ matrix: second sample, observations in rows}
  \item{R}{ number of replicates}
}
\details{
 Computes the coefficient \eqn{\mathcal I}{I_n} and performs a nonparametric
 \eqn{\mathcal E}{E}-test of independence. The test decision is obtained via
 bootstrap, with \code{R} replicates.
 The sample sizes (number of rows) of the two samples must agree, and
 samples must not contain missing values. The statistic
 \eqn{\mathcal E = n \mathcal I^2}{E = I^2} is a ratio of V-statistics based
 on interpoint distances \eqn{\|x_{i}-y_{j}\|}{||x_{i}-y_{j}||}.
 See the reference below for details.
}
\value{
\code{mvI} returns the statistic. \code{mvI.test} returns
 a list with class
 \code{htest} containing
 \item{     method}{ description of test}
 \item{  statistic}{ observed value of the test statistic \eqn{n\mathcal I_n^2}{n I_n^2}}
 \item{   estimate}{ \eqn{\mathcal I_n}{I_n}}
 \item{ replicates}{ replicates of the test statistic}
 \item{    p.value}{ approximate p-value of the test}
 \item{  data.name}{ description of data}
}
\references{
 Bakirov, N.K., Rizzo, M.L., and Szekely, G.J. (2006), A Multivariate
 Nonparametric Test of Independence, \emph{Journal of Multivariate Analysis}
 93/1, 58-80, \cr
 \doi{10.1016/j.jmva.2005.10.005}
 }
 \note{
 Historically this is the first energy test of independence. The
 distance covariance test  \code{\link{dcov.test}}, distance correlation
 \code{\link{dcor}}, and related methods are more recent (2007,2009).
 The distance covariance test is faster and has different properties than
 \code{mvI.test}. Both methods are based on a population independence coefficient
 that characterizes independence and both tests are statistically consistent.
 }
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
 \seealso{
 \code{ \link{indep.test} }
 \code{ \link{mvI.test} }
 \code{ \link{dcov.test} }
 \code{ \link{dcov} }
 }
 \keyword{ htest }
 \keyword{ multivariate }
 \keyword{ nonparametric }
 \concept{ independence }
 \concept{ energy statistics }

