\name{mvI.test}
\alias{mvI.test}
\alias{mvI}
\title{ Independence Coefficient and Test}
\description{
 Computes a type of multivariate nonparametric E-statistic and test of independence
 based on independence coefficient \eqn{\mathcal I_n}{I_n}. This coefficient pre-dates and is different from distance covariance or distance correlation.}
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
 \code{mvI} computes the coefficient \eqn{\mathcal I_n}{I_n} and \code{mvI.test} performs a nonparametric test of independence. The test decision is obtained via permutation 
 bootstrap, with \code{R} replicates.
 The sample sizes (number of rows) of the two samples must agree, and
 samples must not contain missing values. 

Historically this is the first energy test of independence. The
 distance covariance test  \code{\link{dcov.test}}, distance correlation \code{\link{dcor}}, and related methods are more recent (2007, 2009).
 
The distance covariance test \code{\link{dcov.test}} and distance correlation test \code{\link{dcor.test}} are much faster and have  different properties than \code{mvI.test}. All are based on a population independence coefficient that characterizes independence and of these tests are statistically consistent. However, dCor is scale invariant while \eqn{I_n}{I_n} is not. In applications \code{\link{dcor.test}} or \code{\link{dcov.test}} are the recommended tests. 

Computing formula from Bakirov, Rizzo, and Szekely (2006), equation (2): 

Suppose the two samples are \eqn{X_1,\dots,X_n \in R^p} and \eqn{Y_1,\dots,Y_n \in R^q}. Define \eqn{Z_{kl} = (X_k, Y_l) \in R^{p+q}.} 

The independence coefficient \eqn{\mathcal I_n}{I_n} is defined
\deqn{
\mathcal I_n = \sqrt{\frac{2\bar z - z_d - z}{x + y - z}}, 
}
where
\deqn{z_d= \frac{1}{n^2} \sum_{k,l=1}^n |Z_{kk}-Z_{ll}|_{p+q},}
\deqn{z= \frac{1}{n^4} \sum_{k,l=1}^n \sum_{i,j=1}^n |Z_{kl}-Z_{ij}|_{p+q},}
\deqn{\bar z= \frac{1}{n^3} \sum_{k=1}^n \sum_{i,j=1}^n |Z_{kk}-Z_{ij}|_{p+q},}
\deqn{x= \frac{1}{n^2} \sum_{k,l=1}^n |X_{k}-X_{l}|_p,}
\deqn{y= \frac{1}{n^2} \sum_{k,l=1}^n |Y_{k}-Y_{l}|_q.}

Some properties:
\itemize{
\item 
\eqn{0 \leq \mathcal I_n \leq 1} (Theorem 1).
\item 
Large values of \eqn{n \mathcal I_n^2} (or \eqn{\mathcal I_n})  support the alternative hypothesis that the sampled random variables are dependent.
\item \eqn{\mathcal I_n} is invariant to shifts and orthogonal transformations of X and Y.  
\item \eqn{\sqrt{n} \, \mathcal I_n} determines a statistically consistent test of independence against all fixed dependent alternatives (Corollary 1).
\item The population independence coefficient \eqn{\mathcal I} is a normalized distance between the joint characteristic function and the product of the marginal characteristic functions. \eqn{\mathcal I_n} converges almost surely to \eqn{\mathcal I} as \eqn{n \to \infty}. X and Y are independent if and only if \eqn{\mathcal I(X, Y) = 0}.
See the reference below for more details. 
}}
\value{
\code{mvI} returns the statistic. \code{mvI.test} returns
 a list with class
 \code{htest} containing
 \item{     method}{ description of test}
 \item{  statistic}{ observed value of the test statistic \eqn{n\mathcal I_n^2}{n I_n^2}}
 \item{   estimate}{ \eqn{\mathcal I_n}{I_n}}
 \item{ replicates}{ permutation replicates}
 \item{    p.value}{ p-value of the test}
 \item{  data.name}{ description of data}
}
\references{
 Bakirov, N.K., Rizzo, M.L., and Szekely, G.J. (2006), A Multivariate
 Nonparametric Test of Independence, \emph{Journal of Multivariate Analysis} 93/1, 58-80.
 
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.

 Szekely, G.J. and Rizzo, M.L. (2009),
 Brownian Distance Covariance,
 \emph{Annals of Applied Statistics},
 Vol. 3, No. 4, 1236-1265.
 }
\note{
On scale invariance: Distance correlation (\code{\link{dcor}}) has the property that if we change the scale of X from e.g., meters to kilometers, and the scale of Y from e.g. grams to ounces, the statistic and the test are not changed. \eqn{\mathcal I_n}{I_n} does not have this property; it is invariant only under a common rescaling of X and Y by the same constant. Thus, if the units of measurement change for either or both variables, dCor is invariant, but \eqn{\mathcal I_n}{I_n} and possibly the \code{mvI.test}  decision changes.   
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
mvI(iris[1:25, 1], iris[1:25, 2])
\donttest{
mvI.test(iris[1:25, 1], iris[1:25, 2], R=99)
}
}
\seealso{
 \code{ \link{dcov.test} }
 \code{ \link{dcov} }
 \code{ \link{dcor.test} }
 \code{ \link{dcor} }
 \code{ \link{dcov2d} }
 \code{ \link{dcor2d} }
 \code{ \link{indep.test} }
 }
 \keyword{ htest }
 \keyword{ multivariate }
 \keyword{ nonparametric }
 \concept{ independence }
 \concept{ energy statistics }

