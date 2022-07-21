\name{eqdist.etest}
\alias{eqdist.etest}
\alias{eqdist.e}
\alias{ksample.e}
\title{Multisample E-statistic (Energy) Test of Equal Distributions}
\description{
 Performs the nonparametric multisample E-statistic (energy) test
 for equality of multivariate distributions.
}
\usage{
eqdist.etest(x, sizes, distance = FALSE,
    method=c("original","discoB","discoF"), R)
eqdist.e(x, sizes, distance = FALSE,
    method=c("original","discoB","discoF"))
ksample.e(x, sizes, distance = FALSE,
    method=c("original","discoB","discoF"), ix = 1:sum(sizes))
}
\arguments{
  \item{x}{ data matrix of pooled sample}
  \item{sizes}{ vector of sample sizes}
  \item{distance}{logical: if TRUE, first argument is a distance matrix}
  \item{method}{ use original (default) or distance components (discoB, discoF)}
  \item{R}{ number of bootstrap replicates }
  \item{ix}{ a permutation of the row indices of x }
}
\details{
  The k-sample multivariate \eqn{\mathcal{E}}{E}-test of equal distributions
  is performed. The statistic is computed from the original
  pooled samples, stacked in matrix \code{x} where each row
  is a multivariate observation, or the corresponding distance matrix. The
  first \code{sizes[1]} rows of \code{x} are the first sample, the next
  \code{sizes[2]} rows of \code{x} are the second sample, etc.

  The test is implemented by nonparametric bootstrap, an approximate
  permutation test with \code{R} replicates.

  The function \code{eqdist.e} returns the test statistic only; it simply
  passes the arguments through to \code{eqdist.etest} with \code{R = 0}.

    The k-sample multivariate \eqn{\mathcal{E}}{E}-statistic for testing equal distributions
    is returned. The statistic is computed from the original pooled samples, stacked in
    matrix \code{x} where each row is a multivariate observation, or from the distance
    matrix \code{x} of the original data. The
    first \code{sizes[1]} rows of \code{x} are the first sample, the next
    \code{sizes[2]} rows of \code{x} are the second sample, etc.

    The two-sample \eqn{\mathcal{E}}{E}-statistic proposed by
    Szekely and Rizzo (2004)
    is the e-distance \eqn{e(S_i,S_j)}, defined for two samples \eqn{S_i, S_j}
    of size \eqn{n_i, n_j} by
    \deqn{e(S_i,S_j)=\frac{n_i n_j}{n_i+n_j}[2M_{ij}-M_{ii}-M_{jj}],
    }{e(S_i, S_j) = (n_i n_j)(n_i+n_j)[2M_(ij)-M_(ii)-M_(jj)],}
    where
    \deqn{M_{ij}=\frac{1}{n_i n_j}\sum_{p=1}^{n_i} \sum_{q=1}^{n_j}
       \|X_{ip}-X_{jq}\|,}{
       M_{ij} = 1/(n_i n_j) sum[1:n_i, 1:n_j] ||X_(ip) - X_(jq)||,}
       \eqn{\|\cdot\|}{|| ||} denotes Euclidean norm, and \eqn{X_{ip}}{
       X_(ip)} denotes the p-th observation in the i-th sample.


    The original (default method) k-sample
    \eqn{\mathcal{E}}{E}-statistic is defined by summing the pairwise e-distances over
    all \eqn{k(k-1)/2} pairs
    of samples:
    \deqn{\mathcal{E}=\sum_{1 \leq i < j \leq k} e(S_i,S_j).
    }{\emph{E} = sum[i<j] e(S_i,S_j).}
    Large values of \eqn{\mathcal{E}}{\emph{E}} are significant.

    The \code{discoB} method computes the between-sample disco statistic.
    For a one-way analysis, it is related to the original statistic as follows.
    In the above equation, the weights \eqn{\frac{n_i n_j}{n_i+n_j}}{n_i n_j/(n_i+n_j)}
    are replaced with
    \deqn{\frac{n_i + n_j}{2N}\frac{n_i n_j}{n_i+n_j} =
    \frac{n_i n_j}{2N}}{(n_i + n_j)/(2N) n_i n_j/(n_i+n_j) = n_i n_j/(2N)}
    where N is the total number of observations: \eqn{N=n_1+...+n_k}{N=n_1+...+n_k}.

	The \code{discoF} method is based on the disco F ratio, while the \code{discoB}
	method is based on the between sample component.

  Also see \code{disco} and \code{disco.between} functions.
}
\value{
 A list with class \code{htest} containing
 \item{method}{description of test}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{description of data}

 \code{eqdist.e} returns test statistic only.
}
\note{
The pairwise e-distances between samples can be conveniently
computed by the \code{edist} function, which returns a \code{dist} object.
}
\references{
 Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal
 Distributions in High Dimension, \emph{InterStat}, November (5).

 M. L. Rizzo and G. J. Szekely (2010).
 DISCO Analysis: A Nonparametric Extension of
 Analysis of Variance, Annals of Applied Statistics,
 Vol. 4, No. 2, 1034-1055.
\cr \doi{10.1214/09-AOAS245}

 Szekely, G. J. (2000) Technical Report 03-05:
 \eqn{\mathcal{E}}{E}-statistics: Energy of
 Statistical Samples, Department of Mathematics and Statistics, Bowling
 Green State University.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\seealso{
 \code{\link{ksample.e}},
 \code{\link{edist}},
 \code{\link{disco}},
 \code{\link{disco.between}},
 \code{\link{energy.hclust}}.
 }
\examples{
 data(iris)

 ## test if the 3 varieties of iris data (d=4) have equal distributions
 eqdist.etest(iris[,1:4], c(50,50,50), R = 199)

 ## example that uses method="disco"
  x <- matrix(rnorm(100), nrow=20)
  y <- matrix(rnorm(100), nrow=20)
  X <- rbind(x, y)
  d <- dist(X)

  # should match edist default statistic
  set.seed(1234)
  eqdist.etest(d, sizes=c(20, 20), distance=TRUE, R = 199)

  # comparison with edist
  edist(d, sizes=c(20, 10), distance=TRUE)

  # for comparison
  g <- as.factor(rep(1:2, c(20, 20)))
  set.seed(1234)
  disco(d, factors=g, distance=TRUE, R=199)

  # should match statistic in edist method="discoB", above
  set.seed(1234)
  disco.between(d, factors=g, distance=TRUE, R=199)
}
\keyword{ multivariate }
\keyword{ htest }
\keyword{ nonparametric }
\concept{ energy statistics}
