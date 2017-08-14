\name{dcov.test}
\alias{distance covariance}
\alias{dcov.test}
\title{ Distance Covariance Test}
\description{
 Distance covariance test of multivariate independence.
 Distance covariance and distance correlation are
 multivariate measures of dependence.}
\usage{
dcov.test(x, y, index = 1.0, R)
}
\arguments{
  \item{x}{ data or distances of first sample}
  \item{y}{ data or distances of second sample}
  \item{R}{ number of replicates}
  \item{index}{ exponent on Euclidean distance, in (0,2]}
}
\details{
 \code{dcov.test} performs a nonparametric
 test of multivariate independence. The test decision is
 obtained via permutation bootstrap, with \code{R} replicates.

 The sample sizes (number of rows) of the two samples must
 agree, and samples must not contain missing values. Arguments
 \code{x}, \code{y} can optionally be \code{\link{dist}} objects;
 otherwise these arguments are treated as data.

 The statistic is
 \eqn{n \mathcal V_n^2}{nV_n^2} where
 \eqn{\mathcal V_n(x,y)}{V_n(x,y)} = dcov(x,y),
 which is based on interpoint Euclidean distances
 \eqn{\|x_{i}-x_{j}\|}{||x_{i}-x_{j}||}. The \code{index}
 is an optional exponent on Euclidean distance.

Distance correlation is a new measure of dependence between random
vectors introduced by Szekely, Rizzo, and Bakirov (2007).
For all distributions with finite first moments, distance
correlation \eqn{\mathcal R}{R} generalizes the idea of correlation in two
fundamental ways:

 (1) \eqn{\mathcal R(X,Y)}{R(X,Y)} is defined for \eqn{X} and \eqn{Y} in arbitrary dimension.

 (2) \eqn{\mathcal R(X,Y)=0}{R(X,Y)=0} characterizes independence of \eqn{X} and
 \eqn{Y}.

Characterization (2) also holds for powers of Euclidean distance \eqn{\|x_i-x_j\|^s}{|x_i-x_j|^s}, where \eqn{0<s<2}, but (2) does not hold when \eqn{s=2}.

Distance correlation satisfies \eqn{0 \le \mathcal R \le 1}{0 \le R \le 1}, and
\eqn{\mathcal R = 0}{R = 0} only if \eqn{X} and \eqn{Y} are independent. Distance
covariance \eqn{\mathcal V}{V} provides a new approach to the problem of
testing the joint independence of random vectors. The formal
definitions of the population coefficients \eqn{\mathcal V}{V} and
\eqn{\mathcal R}{R} are given in (SRB 2007). The definitions of the
empirical coefficients are given in the energy
\code{\link{dcov}} topic.

For all values of the index in (0,2), under independence
the asymptotic distribution of \eqn{n\mathcal V_n^2}{nV_n^2}
is a quadratic form of centered Gaussian random variables,
with coefficients that depend on the distributions of \eqn{X} and \eqn{Y}. For the general problem of testing independence when the distributions of \eqn{X} and \eqn{Y} are unknown, the test based on \eqn{n\mathcal V^2_n}{n V_n^2} can be implemented as a permutation test. See (SRB 2007) for
theoretical properties of the test, including statistical consistency.
}
\value{
\code{dcov.test} returns a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{   estimate}{ dCov(x,y)}
   \item{  estimates}{ a vector: [dCov(x,y), dCor(x,y), dVar(x), dVar(y)]}
   \item{ replicates}{ replicates of the test statistic}
   \item{    p.value}{ approximate p-value of the test}
   \item{  data.name}{ description of data}
}
\note{For the test of independence,
the distance covariance test statistic is the V-statistic
\eqn{\mathrm{n\, dCov^2} = n \mathcal{V}_n^2}{n V_n^2} (not dCov).
}
\seealso{
 \code{\link{dcov} } \code{\link{dcor} } \code{\link{DCOR}}
 \code{\link{dcor.ttest}}
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

 Szekely, G.J. and Rizzo, M.L. (2009),
 Rejoinder: Brownian Distance Covariance,
 \emph{Annals of Applied Statistics}, Vol. 3, No. 4, 1303-1308.
  }
  \author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- iris[1:50, 1:4]
 y <- iris[51:100, 1:4]
 set.seed(1)
 dcov.test(x, y, R=199)
 \dontrun{
   set.seed(1)
   dcov.test(dist(x), dist(y), R=199)  #same thing
   set.seed(1)
   dcov.test(x, y, index=.5, R=199)
   set.seed(1)
   dcov.test(dist(x), dist(y), index=.5, R=199)  #same thing
 }

 ## Example with dvar=0 (so dcov=0 and pval=1)
 x <- rep.int(1, 10)
 y <- 1:10
 dcov.test(x, y, R=199)
 }
\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

