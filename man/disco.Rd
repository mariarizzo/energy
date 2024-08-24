\name{disco}
\alias{disco}
\alias{disco.between}
\alias{print.disco}
\title{ distance components (DISCO)}
\description{
  E-statistics DIStance COmponents and tests, analogous to variance components
  and anova.
}
\usage{
disco(x, factors, distance, index=1.0, R, method=c("disco","discoB","discoF"))
disco.between(x, factors, distance, index=1.0, R)
}
\arguments{
  \item{x}{ data matrix or distance matrix or dist object}
  \item{factors}{ matrix or data frame of factor labels or integers (not design matrix)}
  \item{distance}{ logical, TRUE if x is distance matrix}
  \item{index}{ exponent on Euclidean distance in (0,2]}
  \item{R}{ number of replicates for a permutation test}
  \item{method}{ test statistic }
}
\details{
  \code{disco} calculates the distance components decomposition of
  total dispersion and if R > 0 tests for significance using the test statistic
  disco "F" ratio (default \code{method="disco"}),
  or using the between component statistic (\code{method="discoB"}),
  each implemented by permutation test.
  
  If \code{x} is a \code{dist} object, argument \code{distance} is
  ignored. If \code{x} is a distance matrix, set \code{distance=TRUE}. 

  In the current release \code{disco} computes the decomposition for one-way models
  only.
}
\value{
  When \code{method="discoF"}, \code{disco} returns a list similar to the
  return value from \code{anova.lm}, and the \code{print.disco} method is
  provided to format the output into a similar table. Details:

\code{disco} returns a class \code{disco} object, which is a list containing
\item{call}{call}
\item{method}{method}
\item{statistic}{vector of observed statistics}
\item{p.value}{vector of p-values}
\item{k}{number of factors}
\item{N}{number of observations}
\item{between}{between-sample distance components}
\item{withins}{one-way within-sample distance components}
\item{within}{within-sample distance component}
\item{total}{total dispersion}
\item{Df.trt}{degrees of freedom for treatments}
\item{Df.e}{degrees of freedom for error}
\item{index}{index (exponent on distance)}
\item{factor.names}{factor names}
\item{factor.levels}{factor levels}
\item{sample.sizes}{sample sizes}
\item{stats}{matrix containing decomposition}

When \code{method="discoB"}, \code{disco} passes the arguments to
\code{disco.between}, which returns a class \code{htest} object.

\code{disco.between} returns a class \code{htest} object, where the test
statistic is the between-sample statistic (proportional to the numerator of the F ratio
of the \code{disco} test.
}
\references{
M. L. Rizzo and G. J. Szekely (2010).
DISCO Analysis: A Nonparametric Extension of
Analysis of Variance, Annals of Applied Statistics,
Vol. 4, No. 2, 1034-1055.
\cr \doi{10.1214/09-AOAS245}
}
\note{
The current version does all calculations via matrix arithmetic and
boot function. Support for more general additive models
and a formula interface is under development.

\code{disco} methods have been added to the cluster distance summary
function \code{edist}, and energy tests for equality of distribution
(see \code{eqdist.etest}).
}
\seealso{
 \code{ \link{edist} }
 \code{ \link{eqdist.e} }
 \code{ \link{eqdist.etest} }
 \code{ \link{ksample.e} }
 }
 \author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and Gabor J. Szekely
 }
\examples{
      ## warpbreaks one-way decompositions
      data(warpbreaks)
      attach(warpbreaks)
      disco(breaks, factors=wool, R=99)
      
      ## warpbreaks two-way wool+tension
      disco(breaks, factors=data.frame(wool, tension), R=0)

      ## warpbreaks two-way wool*tension
      disco(breaks, factors=data.frame(wool, tension, wool:tension), R=0)

      ## When index=2 for univariate data, we get ANOVA decomposition
      disco(breaks, factors=tension, index=2.0, R=99)
      aov(breaks ~ tension)

      ## Multivariate response
      ## Example on producing plastic film from Krzanowski (1998, p. 381)
      tear <- c(6.5, 6.2, 5.8, 6.5, 6.5, 6.9, 7.2, 6.9, 6.1, 6.3,
                6.7, 6.6, 7.2, 7.1, 6.8, 7.1, 7.0, 7.2, 7.5, 7.6)
      gloss <- c(9.5, 9.9, 9.6, 9.6, 9.2, 9.1, 10.0, 9.9, 9.5, 9.4,
                 9.1, 9.3, 8.3, 8.4, 8.5, 9.2, 8.8, 9.7, 10.1, 9.2)
      opacity <- c(4.4, 6.4, 3.0, 4.1, 0.8, 5.7, 2.0, 3.9, 1.9, 5.7,
                   2.8, 4.1, 3.8, 1.6, 3.4, 8.4, 5.2, 6.9, 2.7, 1.9)
      Y <- cbind(tear, gloss, opacity)
      rate <- factor(gl(2,10), labels=c("Low", "High"))

	    ## test for equal distributions by rate
      disco(Y, factors=rate, R=99)
	    disco(Y, factors=rate, R=99, method="discoB")

      ## Just extract the decomposition table
      disco(Y, factors=rate, R=0)$stats

	    ## Compare eqdist.e methods for rate
	    ## disco between stat is half of original when sample sizes equal
	    eqdist.e(Y, sizes=c(10, 10), method="original")
	    eqdist.e(Y, sizes=c(10, 10), method="discoB")

      ## The between-sample distance component
      disco.between(Y, factors=rate, R=0)
}
\keyword{ htest }
\keyword{ multivariate }


