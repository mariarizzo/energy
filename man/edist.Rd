\name{edist}
\alias{edist}
\title{E-distance}
\description{
 Returns the E-distances (energy statistics) between clusters.
}
\usage{
 edist(x, sizes, distance = FALSE, ix = 1:sum(sizes), alpha = 1,
        method = c("cluster","discoB"))
}
\arguments{
  \item{x}{ data matrix of pooled sample or Euclidean distances}
  \item{sizes}{ vector of sample sizes}
  \item{distance}{ logical: if TRUE, x is a distance matrix}
  \item{ix}{ a permutation of the row indices of x }
  \item{alpha}{ distance exponent in (0,2]}
  \item{method}{ how to weight the statistics }
}
\details{
  A vector containing the pairwise two-sample multivariate
  \eqn{\mathcal{E}}{E}-statistics for comparing clusters or samples is returned.
  The e-distance between clusters is computed from the original pooled data,
  stacked in matrix \code{x} where each row is a multivariate observation, or
  from the distance matrix \code{x} of the original data, or distance object
  returned by \code{dist}. The first \code{sizes[1]} rows of the original data
  matrix are the first sample, the next \code{sizes[2]} rows are the second
  sample, etc. The permutation vector \code{ix} may be used to obtain
  e-distances corresponding to a clustering solution at a given level in
  the hierarchy.

  The default method \code{cluster} summarizes the e-distances between
  clusters in a table.
  The e-distance between two clusters \eqn{C_i, C_j}
  of size \eqn{n_i, n_j}
  proposed by Szekely and Rizzo (2005)
  is the e-distance \eqn{e(C_i,C_j)}, defined by
  \deqn{e(C_i,C_j)=\frac{n_i n_j}{n_i+n_j}[2M_{ij}-M_{ii}-M_{jj}],
  }{e(S_i, S_j) = (n_i n_j)/(n_i+n_j)[2M_(ij)-M_(ii)-M_(jj)],}
  where
  \deqn{M_{ij}=\frac{1}{n_i n_j}\sum_{p=1}^{n_i} \sum_{q=1}^{n_j}
     \|X_{ip}-X_{jq}\|^\alpha,}{
     M_{ij} = 1/(n_i n_j) sum[1:n_i, 1:n_j] ||X_(ip) - X_(jq)||^a,}
     \eqn{\|\cdot\|}{|| ||} denotes Euclidean norm, \eqn{\alpha=}{a=}
     \code{alpha}, and \eqn{X_{ip}}{
     X_(ip)} denotes the p-th observation in the i-th cluster.  The
     exponent \code{alpha} should be in the interval (0,2].

  The coefficient \eqn{\frac{n_i n_j}{n_i+n_j}}{(n_i n_j)(n_i+n_j)}
  is one-half of the harmonic mean of the sample sizes. The
  \code{discoB} method is related but with
  different ways of summarizing the pairwise differences between samples.
  The \code{disco} methods apply the coefficient
  \eqn{\frac{n_i n_j}{2N}}{(n_i n_j)/(2N)} where N is the total number
  of observations. This weights each (i,j) statistic by sample size
  relative to N. See the \code{disco} topic for more details.
}
\value{
 A object of class \code{dist} containing the lower triangle of the
 e-distance matrix of cluster distances corresponding to the permutation
 of indices \code{ix} is returned. The \code{method} attribute of the
 distance object is assigned a value of type, index.
}
\references{
 Szekely, G. J. and Rizzo, M. L. (2005) Hierarchical Clustering
 via Joint Between-Within Distances: Extending Ward's Minimum
 Variance Method, \emph{Journal of Classification} 22(2) 151-183.
 \cr \url{http://dx.doi.org/10.1007/s00357-005-0012-9}

M. L. Rizzo and G. J. Szekely (2010).
DISCO Analysis: A Nonparametric Extension of
Analysis of Variance, Annals of Applied Statistics,
Vol. 4, No. 2, 1034-1055.
 \cr \url{http://dx.doi.org/10.1214/09-AOAS245}

 Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal
 Distributions in High Dimension, InterStat, November (5).

 Szekely, G. J. (2000) Technical Report 03-05,
 \eqn{\mathcal{E}}{E}-statistics: Energy of
 Statistical Samples, Department of Mathematics and Statistics,
 Bowling Green State University.
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\seealso{
 \code{\link{energy.hclust}}
 \code{\link{eqdist.etest}}
 \code{\link{ksample.e}}
 \code{\link{disco}}
 }
\examples{
     ## compute cluster e-distances for 3 samples of iris data
     data(iris)
     edist(iris[,1:4], c(50,50,50))
    
     ## pairwise disco statistics
     edist(iris[,1:4], c(50,50,50), method="discoB")  

     ## compute e-distances from a distance object
     data(iris)
     edist(dist(iris[,1:4]), c(50, 50, 50), distance=TRUE, alpha = 1)

     ## compute e-distances from a distance matrix
     data(iris)
     d <- as.matrix(dist(iris[,1:4]))
     edist(d, c(50, 50, 50), distance=TRUE, alpha = 1)

 }
\keyword{ multivariate }
\keyword{ cluster }
\keyword{ nonparametric }
\concept{ energy statistics }
