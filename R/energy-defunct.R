## defunct functions from the energy package

indep.e<-
  function(x, y) {
    # energy statistic for multivariate independence (deprecated)
    .Defunct(new = "mvI", package = "energy")
  }


indep.etest<-
  function(x, y, R) {
    # energy test for multivariate independence (deprecated)
    .Defunct(new = "indep.test", package = "energy",
             msg = "indep.etest removed; use indep.test with method mvI.")
  }


