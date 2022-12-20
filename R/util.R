## util.R
##
## utilities for the energy package
## Author: Maria Rizzo
## github.com/mariarizzo/energy
##


.arg2dist.matrix <- function(x) {
  ## argument check and conversion for energy functions
  ## that take optionally data or distance object arguments
  ## check type of argument, return a distance matrix
  ## supported argument types: matrix, vector, data.frame, tibble, factor, dist
  
  if (anyNA(x)) warning("missing values not supported")
  
  if (inherits(x, "dist")) {
    Dx <- as.matrix(x)
    return(Dx)
  }
  
  if (is.factor(x)) {
    z <- as.matrix(as.integer(x))
    Dx <- calc_dist(z)
    if (!is.ordered(x) && nlevels(x) > 2) {
      # need a 0-1 matrix
      Dx <- matrix(as.integer(Dx > 0), nrow=nrow(Dx))
    }
    return(Dx)
  }
  

  if (is.vector(x) || is.data.frame(x)) {
    ## also for tibble
    Dx <- calc_dist(as.matrix(x))
  }
    
  if (is.matrix(x)) {
    if (is.dmatrix(x)) {
      Dx <- x
    } else {
      ## should be data matrix
      Dx <- calc_dist(x)
    }
  }
  return(Dx)
  
  ## if here, arg type is not supported
  stop(paste("cannot compute distances for", class(x)))
  return(NA)
}

is.dmatrix <- function(x, tol = 100 * .Machine$double.eps) {
  ## check if zero diagonal, symmetric, non-negative square matrix
  ## i.e., distance matrix or dissimilarity matrix
  value <- FALSE
  if (is.matrix(x)) {
    if (nrow(x) == ncol(x)) {
      if (max(abs(diag(x)) < tol) && (max(abs(x - t(x)) < tol))) {
        if (! any(x < 0.0)) value <-  TRUE
      }
    }
  }
  return (value)
}
  
perm.matrix <- function(n, R) {
  ## Generate the same matrix as boot.array with
  ## sim="permutation" and default other arguments
  ## with same seed we get boot.array(boot.out, indices=T)
  
  pfn <- function(x, n) x[sample.int(n)]
  perms <- matrix(1:n, n, R)
  perms <- t(apply(perms, 2, pfn, n=n))
}


permutation <- function(n) {
  ## call the internal permute() function using permute_check()
  J <- 1:n
  a <- .C("permute_check",
          J = as.integer(J),
          n = as.integer(n),
          PACKAGE = "energy")
  return (a$J)
}  

sortrank <- function(x) {
  ## sort and rank data with one call to order()
  ## faster than calling sort and rank separately
  ## returns an object identical to:
  ##   list(x=sort(x), ix=order(x), r=rank(x, ties.method = "first"))
  o <- order(x)
  n <- length(o)
  N <- 1:n
  N[o] <- N
  return(list(x=x[o], ix=o, r=N))
}
