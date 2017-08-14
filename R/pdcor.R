## pdcor.R
##
##

pdcor <- function(x, y, z) {
    if (!(class(x) == "dist")) x <- dist(x)
    if (!(class(y) == "dist")) y <- dist(y)
    if (!(class(z) == "dist")) z <- dist(z)
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    partial_dcor(x, y, z)["pdcor"]
}

pdcov <- function(x, y, z) {
    if (!(class(x) == "dist")) x <- dist(x)
    if (!(class(y) == "dist")) y <- dist(y)
    if (!(class(z) == "dist")) z <- dist(z)
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    partial_dcov(x, y, z)
}

