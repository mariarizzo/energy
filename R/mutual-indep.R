mutualIndep.test <-
function(x, R) {
  if (NCOL(x) < 2) {
    stop("Expecting two or more samples")
    }
  bootfn <- function(x, i) { 
    d <- ncol(x)
    dc <- numeric(d-1)
    for (k in 1:(d-1)) {
      dc[k] <- energy::bcdcor(x[i,k], x[,(k+1):d])
    }
    return (dc)
  }
  
  b <- boot::boot(x, bootfn, sim="permutation", R=R)
  t0 <- sum(b$t0)
  tp <- rowSums(b$t)
  pval <- (1 + sum(tp > t0)) / (R + 1)
  estimate <- round(b$t0, 3)
  names(t0) <- "Sum(R*)"
  names(estimate) <- paste0("R*", 1:length(b$t0))
  method <- paste("Energy Test of Mutual Independence")
  call <- match.call()
  NOTE <- "statistic=sum(bcdcor); permutation test"
  rval <- list(statistic = t0, p.value = pval, call = call,
               data.name=paste(deparse(substitute(x))," dim ", paste(dim(x), collapse=",")),
               estimate=estimate, method=method, note=NOTE)
  class(rval) <- "power.htest"
  return(rval)
}

