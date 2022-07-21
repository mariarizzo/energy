\name{dcov.test}
\alias{distance covariance}
\alias{dcov.test}
\alias{dcor.test}
\title{ Distance Covariance Test and Distance Correlation test}
\description{
 Distance covariance test and distance correlation test of multivariate independence.
 Distance covariance and distance correlation are
 multivariate measures of dependence.}
\usage{
dcov.test(x, y, index = 1.0, R = NULL)
dcor.test(x, y, index = 1.0, R)
}
\arguments{
  \item{x}{ data or distances of first sample}
  \item{y}{ data or distances of second sample}
  \item{R}{ number of replicates}
  \item{index}{ exponent on Euclidean distance, in (0,2]}
}
\details{
 \code{dcov.test} and \code{dcor.test} are nonparametric
 tests of multivariate independence. The test decision is
 obtained via permutation bootstrap, with \code{R} replicates.

 The sample sizes (number of rows) of the two samples must
 agree, and samples must not contain missing values. 
 
 The \code{index} is an optional exponent on Euclidean distance.
Valid exponents for energy are in (0, 2) excluding 2. 

Argument types supported are 
numeric data matrix, data.frame, or tibble, with observations in rows;
numeric vector; ordered or unordered factors. In case of unordered factors
a 0-1 distance matrix is computed.

Optionally pre-computed distances can be input as class "dist" objects or as distance matrices. 
For data types of arguments,
distance matrices are computed internally. 

 The \code{dcov} test statistic is
 \eqn{n \mathcal V_n^2}{nV_n^2} where
 \eqn{\mathcal V_n(x,y)}{V_n(x,y)} = dcov(x,y),
 which is based on interpoint Euclidean distances
 \eqn{\|x_{i}-x_{j}\|}{||x_{i}-x_{j}||}. The \code{index}
 is an optional exponent on Euclidean distance.

Similarly, the \code{dcor} test statistic is based on the normalized
coefficient, the distance correlation. (See the manual page for \code{dcor}.)

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
\code{dcov.test} or \code{dcor.test} returns a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{   estimate}{ dCov(x,y) or dCor(x,y)}
   \item{  estimates}{ a vector: [dCov(x,y), dCor(x,y), dVar(x), dVar(y)]}
   \item{  condition}{ logical, permutation test applied}
   \item{ replicates}{ replicates of the test statistic}
   \item{    p.value}{ approximate p-value of the test}
   \item{          n}{ sample size}
   \item{  data.name}{ description of data}
}
\note{For the dcov test of independence,
the distance covariance test statistic is the V-statistic
\eqn{\mathrm{n\, dCov^2} = n \mathcal{V}_n^2}{n V_n^2} (not dCov). 
}
\seealso{
 \code{\link{dcov} } \code{\link{dcor} } 
 \code{\link{pdcov.test}} \code{\link{pdcor.test}} 
 \code{\link{dcor.ttest}}
 }

\references{
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \doi{10.1214/009053607000000505}

 Szekely, G.J. and Rizzo, M.L. (2009),
 Brownian Distance Covariance,
 \emph{Annals of Applied Statistics},
 Vol. 3, No. 4, 1236-1265.
 \cr \doi{10.1214/09-AOAS312}

 Szekely, G.J. and Rizzo, M.L. (2009),
 Rejoinder: Brownian Distance Covariance,
 \emph{Annals of Applied Statistics}, Vol. 3, No. 4, 1303-1308.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- iris[1:50, 1:4]
 y <- iris[51:100, 1:4]
 set.seed(1)
 dcor.test(dist(x), dist(y), R=199)
 set.seed(1)
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

