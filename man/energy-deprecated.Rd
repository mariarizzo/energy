\name{energy-deprecated}
\alias{dcor.ttest}
\alias{dcor.t}
\title{ Deprecated Functions}
\description{ These deprecated functions have been replaced by revised functions and will be removed in future releases of the energy package.}
\usage{
dcor.ttest(x, y, distance=FALSE)
dcor.t(x, y, distance=FALSE)
}
\arguments{
  \item{x}{ data or distances of first sample}
  \item{y}{ data or distances of second sample}
  \item{distance}{ logical: TRUE if x and y are distances}
}
\details{
  * dcor.t has been replaced by \code{dcorT}. See \link{dcorT} for details.
  * dcor.ttest has been replaced by \code{dcorT.test}. See \link{dcorT.test} for details. 
}
\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

