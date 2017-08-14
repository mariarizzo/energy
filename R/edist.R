edist <-
function(x, sizes, distance = FALSE, ix = 1:sum(sizes), alpha = 1,
    method = c("cluster","discoB","discoF")) {
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
    if (distance == FALSE) {
        if (is.vector(x)) x <- matrix(x, nrow = length(x), ncol = 1)
        dst <- as.matrix(dist(x))
        }
    else dst <- as.matrix(x)
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


    if (type == "discoF" || type == "discoB") {
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
            N <- NROW(d)
	    total <- sum(d) / (2*N)
	    trt <- factor(c(rep(1,n1),rep(2,n2)))
	    y <- as.vector(d[,1])
	    M <- model.matrix(y ~ 0 + trt)
	    G <- t(M) %*% d %*% M
	    withins <- diag(G) / (2*c(n1,n2))
	    W <- sum(withins)
	    B <- total - W
	    ifelse (type == "discoF",
		e[i,j] <- e[j,i] <- B / (W/(N-2)),
		e[i,j] <- e[j,i] <- B)
            }
        }
	}
    e <- as.dist(e)
    attr(e,"method") <- paste(method,": index= ", alpha)
    e
}

