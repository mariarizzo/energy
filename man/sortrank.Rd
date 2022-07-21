\name{sortrank}
\alias{sortrank}
\title{ Sort, order and rank a vector }
\description{
A utility that returns a list with the components
equivalent to sort(x), order(x), rank(x, ties.method = "first").
 }
\usage{
sortrank(x)
}
\arguments{
  \item{x}{ vector compatible with sort(x)}
}
\details{
This utility exists to save a little time on large vectors when two or all three of the sort(), order(), rank() results are required. In case of ties, the ranks component matches \code{rank(x, ties.method = "first")}. 
}
\value{
A list with components 
  \item{x}{the sorted input vector x}
  \item{ix}{the permutation = order(x) which rearranges x into ascending order}
  \item{r}{the ranks of x}
}
\note{
This function was benchmarked faster than the combined calls to \code{sort} and \code{rank}. 
}
\examples{
sortrank(rnorm(5))
}
\references{
See \code{\link{sort}}.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu}
}

