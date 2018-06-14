## util.R
##
## miscellaneous utilities
##


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

