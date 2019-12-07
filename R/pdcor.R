## pdcor.R
##
##

pdcor <- function(x, y, z) {
    if (!inherits(x, "dist")) x <- dist(x)
    if (!inherits(y, "dist")) y <- dist(y)
    if (!inherits(z, "dist")) z <- dist(z)
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    partial_dcor(x, y, z)["pdcor"]
}

pdcov <- function(x, y, z) {
    if (!inherits(x, "dist")) x <- dist(x)
    if (!inherits(y, "dist")) y <- dist(y)
    if (!inherits(z, "dist")) z <- dist(z)
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    partial_dcov(x, y, z)
}

