\name{EVnormal}
\docType{data}
\alias{EVnormal}
\alias{eigenvalues}
\title{Eigenvalues for the energy Test of Univariate Normality}
\description{
  Pre-computed eigenvalues corresponding to the asymptotic sampling
  distribution of the energy test statistic for univariate
  normality, under the null hypothesis. Four Cases are computed:
    \enumerate{
      \item Simple hypothesis, known parameters.
      \item Estimated mean, known variance.
      \item Known mean, estimated variance.
      \item Composite hypothesis, estimated parameters.
    }
Case 4 eigenvalues are used in the test function \code{normal.test}
when \code{method=="limit"}.
}
\usage{data(EVnormal)}
\format{Numeric matrix with 125 rows and 5 columns;
  column 1 is the index, and columns 2-5 are 
  the eigenvalues of Cases 1-4.}
\source{Computed}
\references{
  Szekely, G. J. and Rizzo, M. L. (2005) A New Test for
  Multivariate Normality, \emph{Journal of Multivariate Analysis},
  93/1, 58-80,
  \doi{10.1016/j.jmva.2003.12.002}.
}
