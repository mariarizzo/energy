# deprecated independence test

indep.test<-
function(x, y, method = c("dcov","mvI"), index = 1, R) {
  # two energy tests for multivariate independence
  .Deprecated(new = "dcov.test", package = "energy",
              msg = "indep.test is deprecated, 
              replaced by dcov.test or mvI.test")
  
  type <- match.arg(method)
    if (type == "dcov")
        return(dcov.test(x, y, index, R)) else
    if (type == "mvI")
        return(mvI.test(x, y, R))
}
