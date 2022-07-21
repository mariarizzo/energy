\name{distance correlation}
\alias{dcor}
\alias{dcov}
\title{ Distance Correlation and Covariance Statistics}
\description{
 Computes distance covariance and distance correlation statistics,
 which are multivariate measures of dependence.
 }
\usage{
dcov(x, y, index = 1.0)
dcor(x, y, index = 1.0)
}
\arguments{
  \item{x}{ data or distances of first sample}
  \item{y}{ data or distances of second sample}
  \item{index}{ exponent on Euclidean distance, in (0,2]}
}
\details{
 \code{dcov} and \code{dcor} compute distance
 covariance and distance correlation statistics.

  The sample sizes (number of rows) of the two samples must
  agree, and samples must not contain missing values. 
  
The \code{index} is an optional exponent on Euclidean distance.
Valid exponents for energy are in (0, 2) excluding 2. 

Argument types supported are 
numeric data matrix, data.frame, or tibble, with observations in rows;
numeric vector; ordered or unordered factors. In case of unordered factors
a 0-1 distance matrix is computed.

Optionally pre-computed distances can be input as class "dist" objects or as distance matrices. 
 For data types of arguments, distance matrices are computed internally. 

Distance correlation is a new measure of dependence between random
vectors introduced by Szekely, Rizzo, and Bakirov (2007).
For all distributions with finite first moments, distance
correlation \eqn{\mathcal R}{R} generalizes the idea of correlation in two
fundamental ways:
 (1) \eqn{\mathcal R(X,Y)}{R(X,Y)} is defined for \eqn{X} and \eqn{Y} in arbitrary dimension.
 (2) \eqn{\mathcal R(X,Y)=0}{R(X,Y)=0} characterizes independence of \eqn{X} and
 \eqn{Y}.

Distance correlation satisfies \eqn{0 \le \mathcal R \le 1}{0 \le R \le 1}, and
\eqn{\mathcal R = 0}{R = 0} only if \eqn{X} and \eqn{Y} are independent. Distance
covariance \eqn{\mathcal V}{V} provides a new approach to the problem of
testing the joint independence of random vectors. The formal
definitions of the population coefficients \eqn{\mathcal V}{V} and
\eqn{\mathcal R}{R} are given in (SRB 2007). The definitions of the
empirical coefficients are as follows.

The empirical distance covariance \eqn{\mathcal{V}_n(\mathbf{X,Y})}{V_n(X,Y)}
with index 1 is
the nonnegative number defined by
\deqn{
 \mathcal{V}^2_n (\mathbf{X,Y}) = \frac{1}{n^2} \sum_{k,\,l=1}^n
 A_{kl}B_{kl}
 }{
V^2_n (X,Y) = (1/n^2) sum_{k,l=1:n}
 A_{kl}B_{kl}
 }
 where \eqn{A_{kl}} and \eqn{B_{kl}} are
 \deqn{
A_{kl} = a_{kl}-\bar a_{k.}- \bar a_{.l} + \bar a_{..}
}
\deqn{
 B_{kl} = b_{kl}-\bar b_{k.}- \bar b_{.l} + \bar b_{..}.
 }
Here
\deqn{
a_{kl} = \|X_k - X_l\|_p, \quad b_{kl} = \|Y_k - Y_l\|_q, \quad
k,l=1,\dots,n,
}{
a_{kl} = ||X_k - X_l||_p, b_{kl} = ||Y_k - Y_l||_q,
k,l=1,\dots,n,
}
and the subscript \code{.} denotes that the mean is computed for the
index that it replaces.  Similarly,
\eqn{\mathcal{V}_n(\mathbf{X})}{V_n(X)} is the nonnegative number defined by
\deqn{
 \mathcal{V}^2_n (\mathbf{X}) = \mathcal{V}^2_n (\mathbf{X,X}) =
 \frac{1}{n^2} \sum_{k,\,l=1}^n
 A_{kl}^2.
 }{
V^2_n (X) = V^2_n (X,X) =
 (1/n^2) sum_{k,l=1:n}
 A_{kl}^2.
 }

The empirical distance correlation \eqn{\mathcal{R}_n(\mathbf{X,Y})}{R(\mathbf{X,Y})} is
the square root of
\deqn{
  \mathcal{R}^2_n(\mathbf{X,Y})=
 \frac {\mathcal{V}^2_n(\mathbf{X,Y})}
 {\sqrt{ \mathcal{V}^2_n (\mathbf{X}) \mathcal{V}^2_n(\mathbf{Y})}}.
}{
  R^2_n(X,Y)=
 V^2_n(X,Y) / sqrt(V^2_n (X) V^2_n(Y)).
}
See \code{\link{dcov.test}} for a test of multivariate independence
based on the distance covariance statistic. 
}
\value{
\code{dcov} returns the sample distance covariance and
\code{dcor} returns the sample distance correlation.
}
\note{
Note that it is inefficient to compute dCor by:

square root of
\code{dcov(x,y)/sqrt(dcov(x,x)*dcov(y,y))}

because the individual
calls to \code{dcov} involve unnecessary repetition of calculations.
}
\seealso{
\code{\link{dcov2d}} \code{\link{dcor2d}} 
\code{\link{bcdcor}}  \code{\link{dcovU}}  \code{\link{pdcor}}
\code{\link{dcov.test}} \code{\link{dcor.test}}  \code{\link{pdcor.test}}
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
 dcov(x, y)
 dcov(dist(x), dist(y))  #same thing
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

