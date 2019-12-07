
kgroups <- function(x, k, iter.max = 10, nstart = 1, cluster = NULL) {
  distance <- inherits(x, "dist")
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("x must be numeric")
  n <- nrow(x)
  if (is.null(cluster)) {
    cluster <- sample(0:(k-1), size = n, replace = TRUE)
  } else {
    ## recode cluster as 0,1,...,k-1
    cluster <- factor(cluster)
    if(length(levels(cluster)) != k)
      stop("cluster vector does not have k clusters")
    cluster <- as.integer(cluster) - 1
    if(length(cluster) != n)
      stop("data and length of cluster vector must match")
  }
  value <- kgroups_start(x, k, cluster, iter.max, distance = distance)

  if (nstart > 1) {
    objective <- rep(0, nstart)
    objective[1] <- value$W
    values <- vector("list", nstart)
    values[[1]] <- value
    for (j in 2:nstart) {
      ## random initialization of cluster labels
      cluster <- sample(0:(k-1), size = n, replace = TRUE)
      values[[j]] <- kgroups_start(x, k, cluster, iter.max, distance = distance)

      objective[j] <- values[[j]]$W
    }
    best <- which.min(objective)
    value <- values[[best]]
  }

  obj  <- structure(list(
    call = match.call(),
    cluster = value$cluster + 1,
    sizes = value$sizes,
    within = value$within,
    W = sum(value$within),
    count = value$count,
    iterations = value$it,
    k = k),
    class = "kgroups")
  return (obj)
}


print.kgroups <- function(x, ...) {
  cat("\n"); print(x$call)
  cat("\nK-groups cluster analysis\n")
  cat(x$k, " groups of size ", x$sizes, "\n")
  cat("Within cluster distances:\n", x$within)
  cat("\nIterations: ", x$iterations, "  Count: ", x$count, "\n")
}

fitted.kgroups <- function(object, method = c("labels", "groups"), ...) {
  method = match.arg(method)
  if (method == "groups") {
    k <- object$k
    CList <- vector("list", k)
    for (i in 1:k)
      CList[[i]] <- which(object$cluster == i)
    return (CList)
  }
  return (object$cluster)
}
