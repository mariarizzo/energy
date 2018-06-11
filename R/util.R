## util.R
##
## miscellaneous utilities
##


sortrank <- function(x) {
  ## sort and rank data with one call to order()
  ## faster than sort(x, index.return=TRUE) because
  ## x[order(x)] is faster than sort(x)
  ## returns an object identical to:
  ##   list(x=sort(x), ix=order(x), r=rank(x))
  o <- order(x)
  n <- length(o)
  N <- 1:n
  N[o] <- N
  return(list(x=x[o], ix=o, r=N))
}

