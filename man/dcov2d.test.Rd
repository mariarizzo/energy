\name{dcov2d.test}
\alias{dcov2d.test}
\title{ Large-Sample Distance Covariance Test}
\description{
A large-sample distance covariance test of independence is implemented for the O(n log n) V-statistic \code{dcov2d}. The computing algorithm for \code{dcov2d} is implemented for real valued X and Y (bivariate data).
}
\usage{
dcov2d.test(x, y, m=200)
}
\arguments{
  \item{x}{ data vector of first sample}
  \item{y}{ data vector of second sample}
  \item{m}{ number of eigenvalues}
}
\details{
For details about testing independence with distance covariance see \code{dcov}, 
 \code{dcov.test} and \code{dcor}. 
 
 \code{dcov2d.test} implements a large-sample test of independence for the V-statistic returned by \code{dcov2d}, which is based on a fast O(n log n) computing algorithm. The test applies the limit distribution of the test statistic \eqn{n V_n}. The asymptotic distribution of the test statistic under independence is of the form 
 \deqn{
 \sum_{i=1}^{\infty} \lambda_i Z_i^2,
 }{\sum_{i=1}^{\infty} \lambda_i Z_i^2,}
 where \eqn{\lambda_i > 0}{\lambda_i > 0} are eigenvalues and \eqn{Z_i} are iid standard normal random variables.  
 The eigenvalues are estimated from a subsample of the data, and the upper tail probabilities are computed by Imhof's method. 
 
The argument \code{m} specifies the number of eigenvalues estimated for computing the asymptotic null distribution of the test statistic. Use the permutation test \code{dcov.test} or \code{dcor.test} if n < 200.   
}
\value{
\code{dcov2d.test} or \code{dcor.test} returns a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{   estimate}{ V_n, R_n (squared dcor) }
   \item{    p.value}{ approximate p-value of the test}
   \item{     lambda}{ eigenvalues}
   \item{  data.name}{ description of data}
}
\seealso{
 \code{\link{dcov2d}} \code{\link{dcor2d}} (fast computing formulas)
 
 \code{\link{dcov}} \code{\link{dcov.test}}  \code{\link{dcor}} \code{\link{dcor.test}} (multivariate statistics and permutation test) 
}

\references{
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \url{http://dx.doi.org/10.1214/009053607000000505}

 Szekely, G.J. and Rizzo, M.L. (2009),
 Brownian Distance Covariance,
 \emph{Annals of Applied Statistics},
 Vol. 3, No. 4, 1236-1265.
 \cr \url{http://dx.doi.org/10.1214/09-AOAS312}
 
Huo, X. and Szekely, G.J. (2016). Fast computing for 
distance covariance. Technometrics, 58(4), 435-447.
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
   \donttest{
    x <- rlnorm(400)
    y <- rexp(400)
    dcov.test(x, y, R=199)
    dcov2d.test(x, y)
    }
}
\keyword{ htest }
\keyword{ nonparametric }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

