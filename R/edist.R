edist <-
function(x, sizes, distance = FALSE, ix = 1:sum(sizes), alpha = 1,
    method = c("cluster","discoB")) {
    #  computes the e-dissimilarity matrix between k samples or clusters
    #  x:          pooled sample or Euclidean distances
    #  sizes:      vector of sample (cluster) sizes
    #  distance:   TRUE if x is a distance matrix, otherwise FALSE
    #  ix:         a permutation of row indices of x
    #  alpha:      distance exponent
    #  method:     cluster distances or disco statistics
    #
    k <- length(sizes)
    if (k == 1) return (as.dist(0.0))
    if (k < 1) return (NA)
    e <- matrix(nrow=k, ncol=k)
    n <- cumsum(sizes)
    m <- 1 + c(0, n[1:(k-1)])
    
    if (is.vector(x)) x <- matrix(x, ncol=1)
    if (inherits(x, "dist")) distance <- TRUE
    if (distance)
      dst <- as.matrix(x) else dst <- as.matrix(dist(x))
    N <- NROW(dst)
    if (NCOL(dst) != N)
      stop("distance==TRUE but first argument is not distance")
    
    if (alpha != 1) {
    	if (alpha <= 0 || alpha > 2)
    	    warning("exponent alpha should be in (0,2]")
    	dst <- dst^alpha
    	}

    type <- match.arg(method)
	  if (type == "cluster") {
      for (i in 1:(k - 1)) {
        e[i, i] <- 0.0
        for (j in (i + 1):k) {
            n1 <- sizes[i]
            n2 <- sizes[j]
            ii <- ix[m[i]:n[i]]
            jj <- ix[m[j]:n[j]]
            w <- n1 * n2 / (n1 + n2)
            m11 <- sum(dst[ii, ii]) / (n1 * n1)
            m22 <- sum(dst[jj, jj]) / (n2 * n2)
            m12 <- sum(dst[ii, jj]) / (n1 * n2)
            e[i, j] <- e[j, i] <- w * ((m12 + m12) - (m11 + m22))
            }
        }
    }


    if (type == "discoB") {
      #disco statistics for testing F=G
      for (i in 1:(k - 1)) {
        e[i, i] <- 0.0
        for (j in (i + 1):k) {
            n1 <- sizes[i]
            n2 <- sizes[j]
            ii <- ix[m[i]:n[i]]
            jj <- ix[m[j]:n[j]]
            J <- c(ii,jj)
            d <- dst[J, J]
      	    e[i, j] <- eqdist.e(d, sizes=c(n1, n2), distance=TRUE)
      	    e[j, i] <- e[i, j] <- e[i, j] * (n1 + n2) 
        }
      }
      e <- 0.5 * e / sum(sizes)  #discoB formula
    }
    
    e <- as.dist(e) 
    attr(e,"method") <- paste(method,": index= ", alpha)
    e
}

