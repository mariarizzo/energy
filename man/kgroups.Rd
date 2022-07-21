\name{kgroups}
\alias{kgroups}

\title{
K-Groups Clustering
}
\description{
Perform k-groups clustering by energy distance.
}
\usage{
kgroups(x, k, iter.max = 10, nstart = 1, cluster = NULL)
}
\arguments{
  \item{x}{Data frame or data matrix or distance object}
  \item{k}{number of clusters}
  \item{iter.max}{maximum number of iterations}
  \item{nstart}{number of restarts}
  \item{cluster}{initial clustering vector}
}

\details{
K-groups is based on the multisample energy distance for comparing distributions.
Based on the disco decomposition of total dispersion (a Gini type mean distance) the objective function should either maximize the total between cluster energy distance, or equivalently, minimize the total within cluster energy distance. It is more computationally efficient to minimize within distances, and that makes it possible to use a modified version of the Hartigan-Wong algorithm (1979) to implement K-groups clustering.

The within cluster Gini mean distance is
\deqn{G(C_j) = \frac{1}{n_j^2} \sum_{i,m=1}^{n_j} |x_{i,j} - x_{m,j}|}
and the K-groups within cluster distance is
\deqn{W_j = \frac{n_j}{2}G(C_j) = \frac{1}{2 n_j} \sum_{i,m=1}^{n_j} |x_{i,j} - x_{m,j}.}
If z is the data matrix for cluster \eqn{C_j}, then \eqn{W_j} could be computed as
\code{sum(dist(z)) / nrow(z)}.

If cluster is not NULL, the clusters are initialized by this vector (can be a factor or integer vector). Otherwise clusters are initialized with random labels in k approximately equal size clusters.

If \code{x} is not a distance object (class(x) == "dist") then \code{x} is converted to a data matrix for analysis. 

Run up to \code{iter.max} complete passes through the data set until a local min is reached. If \code{nstart > 1}, on second and later starts, clusters are initialized at random, and the best result is returned. 
}

\value{
An object of class \code{kgroups} containing the components
\item{call}{the function call}
\item{cluster}{vector of cluster indices}
\item{sizes}{cluster sizes}
\item{within}{vector of Gini within cluster distances}
\item{W}{sum of within cluster distances}
\item{count}{number of moves}
\item{iterations}{number of iterations}
\item{k}{number of clusters}

\code{cluster} is a vector containing the group labels, 1 to k. \code{print.kgroups}
prints some of the components of the kgroups object.

Expect that count is 0 if the algorithm converged to a local min (that is, 0 moves happened on the last iteration). If iterations equals iter.max and count is positive, then the algorithm did not converge to a local min. 
}
\author{
Maria Rizzo and Songzi Li
}

\references{
Li, Songzi (2015).
"K-groups: A Generalization of K-means by Energy Distance."
Ph.D. thesis, Bowling Green State University.

Li, S. and Rizzo, M. L. (2017).
"K-groups: A Generalization of K-means Clustering".
ArXiv e-print 1711.04359. https://arxiv.org/abs/1711.04359

Szekely, G. J., and M. L. Rizzo. "Testing for equal distributions in high dimension." InterStat 5, no. 16.10 (2004).

Rizzo, M. L., and G. J. Szekely. "Disco analysis: A nonparametric extension of analysis of variance." The Annals of Applied Statistics (2010): 1034-1055.

Hartigan, J. A. and Wong, M. A. (1979). "Algorithm AS 136: A K-means clustering algorithm." Applied Statistics, 28, 100-108. doi: 10.2307/2346830.
}

\examples{
  x <- as.matrix(iris[ ,1:4])
  set.seed(123)
  kg <- kgroups(x, k = 3, iter.max = 5, nstart = 2)
  kg
  fitted(kg)
  
  \donttest{
    d <- dist(x)
    set.seed(123)
    kg <- kgroups(d, k = 3, iter.max = 5, nstart = 2)
    kg
    
    kg$cluster
  
    fitted(kg)
    fitted(kg, method = "groups")
    }
}

\keyword{ cluster }
\keyword{ multivariate }
