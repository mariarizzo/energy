\name{energy.hclust}
\alias{energy.hclust}
\title{ Hierarchical Clustering by Minimum (Energy) E-distance }
\description{
  Performs hierarchical clustering by minimum (energy) E-distance method.
}
\usage{
    energy.hclust(dst, alpha = 1)
}
\arguments{
  \item{dst}{\code{dist} object}
  \item{alpha}{distance exponent}
}
\details{
  Dissimilarities are \eqn{d(x,y) = \|x-y\|^\alpha}{||x-y||^a},
  where the exponent \eqn{\alpha}{a} is in the interval (0,2].
  This function performs agglomerative hierarchical clustering.
  Initially, each of the n singletons is a cluster. At each of n-1 steps, the
  procedure merges the pair of clusters with minimum e-distance.
  The e-distance between two clusters \eqn{C_i, C_j} of sizes \eqn{n_i, n_j} 
  is given by
    \deqn{e(C_i, C_j)=\frac{n_i n_j}{n_i+n_j}[2M_{ij}-M_{ii}-M_{jj}],
    }
    where
    \deqn{M_{ij}=\frac{1}{n_i n_j}\sum_{p=1}^{n_i} \sum_{q=1}^{n_j}
       \|X_{ip}-X_{jq}\|^\alpha,}{
       M_{ij} = 1/(n_i n_j) sum[1:n_i, 1:n_j] ||X_(ip) - X_(jq)||^a,}
       \eqn{\|\cdot\|}{|| ||} denotes Euclidean norm, and \eqn{X_{ip}}{
       X_(ip)} denotes the p-th observation in the i-th cluster.

  The return value is an object of class \code{hclust}, so \code{hclust}
  methods such as print or plot methods, \code{plclust}, and \code{cutree}
  are available. See the documentation for \code{hclust}.

  The e-distance measures both the heterogeneity between clusters and the
  homogeneity within clusters. \eqn{\mathcal E}{E}-clustering
  (\eqn{\alpha=1}{a=1}) is particularly effective in
  high dimension, and is more effective than some standard hierarchical
  methods when clusters have equal means (see example below).
  For other advantages see the references.
  
  \code{edist} computes the energy distances for the result (or any partition)
  and returns the cluster distances in a \code{dist} object. See the \code{edist}
  examples.
}
\value{
     An object of class \code{hclust} which describes the tree produced by
     the clustering process. The object is a list with components:
\item{merge:}{ an n-1 by 2 matrix, where row i of \code{merge} describes the
     merging of clusters at step i of the clustering. If an element j in the
     row is negative, then observation -j was merged at this
     stage. If j is positive then the merge was with the cluster
     formed at the (earlier) stage j of the algorithm.}
\item{height:}{the clustering height: a vector of n-1 non-decreasing
     real numbers (the e-distance between merging clusters)}
\item{order:}{ a vector giving a permutation of the indices of
     original observations suitable for plotting, in the sense that a
     cluster plot using this ordering and matrix \code{merge} will not have
     crossings of the branches.}
\item{labels:}{ labels for each of the objects being clustered.}
\item{call:}{ the call which produced the result.}
\item{method:}{ the cluster method that has been used (e-distance).}
\item{dist.method:}{ the distance that has been used to create \code{dst}.}
}
\note{
Currently \code{stats::hclust} implements Ward's method by \code{method="ward.D2"},
which applies the squared distances. That method was previously \code{"ward"}. 
Because both \code{hclust} and energy use the same type of Lance-Williams recursive formula to update cluster distances, now with the additional option \code{method="ward.D"} in \code{hclust}, the
energy distance method is easily implemented by \code{hclust}. (Some "Ward" algorithms do not use Lance-Williams, however). Energy clustering (with \code{alpha=1}) and "ward.D" now return the same result, except that the cluster heights of energy hierarchical clustering with \code{alpha=1} are two times the heights from \code{hclust}. In order to ensure compatibility with hclust methods, \code{energy.hclust} now passes arguments through to \code{hclust} after possibly applying the optional exponent to distance.
}
\references{
     Szekely, G. J. and Rizzo, M. L. (2005) Hierarchical Clustering
     via Joint Between-Within Distances: Extending Ward's Minimum
     Variance Method, \emph{Journal of Classification} 22(2) 151-183.
     \cr \doi{10.1007/s00357-005-0012-9}

     Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal
     Distributions in High Dimension, \emph{InterStat}, November (5).

 Szekely, G. J. (2000) Technical Report 03-05:
 \eqn{\mathcal{E}}{E}-statistics: Energy of
 Statistical Samples, Department of Mathematics and Statistics, Bowling
 Green State University.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
  Gabor J. Szekely
  }
\seealso{
  \code{\link{edist}} \code{\link{ksample.e}} \code{\link{eqdist.etest}} \code{hclust}}
\examples{
   \dontrun{
   library(cluster)
   data(animals)
   plot(energy.hclust(dist(animals)))

   data(USArrests)
   ecl <- energy.hclust(dist(USArrests))
   print(ecl)
   plot(ecl)
   cutree(ecl, k=3)
   cutree(ecl, h=150)

   ## compare performance of e-clustering, Ward's method, group average method
   ## when sampled populations have equal means: n=200, d=5, two groups
   z <- rbind(matrix(rnorm(1000), nrow=200), matrix(rnorm(1000, 0, 5), nrow=200))
   g <- c(rep(1, 200), rep(2, 200))
   d <- dist(z)
   e <- energy.hclust(d)
   a <- hclust(d, method="average")
   w <- hclust(d^2, method="ward.D2")
   list("E" = table(cutree(e, k=2) == g), "Ward" = table(cutree(w, k=2) == g),
    "Avg" = table(cutree(a, k=2) == g))
  }
 }
\keyword{ multivariate }
\keyword{ cluster }
\concept{ energy statistics }
