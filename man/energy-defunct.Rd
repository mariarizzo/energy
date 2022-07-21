\name{indep.etest}
\alias{indep.e}
\alias{indep.etest}
\title{ Energy Statistic Test of Independence}
\description{Defunct: use \code{indep.test} with \code{method = mvI}.
 Computes a multivariate nonparametric E-statistic and test of independence.}
\usage{
indep.e(x, y)
indep.etest(x, y, R)
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
 The sample coefficient \eqn{\mathcal I}{I} is returned by \code{indep.e}.
 The function \code{indep.etest} returns a list with class
 \code{htest} containing
 \item{method}{description of test}
 \item{statistic}{observed value of the coefficient \eqn{\mathcal I}{I}}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{description of data}
}
\references{
 Bakirov, N.K., Rizzo, M.L., and Szekely, G.J. (2006), A Multivariate
 Nonparametric Test of Independence, \emph{Journal of Multivariate Analysis}
 93/1, 58-80
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}

\keyword{ htest }
\keyword{ multivariate }
\concept{ energy statistics }

