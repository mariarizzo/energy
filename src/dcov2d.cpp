#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.calcH2d)]]
NumericMatrix calcH2d(NumericVector x, NumericVector y, IntegerVector I,
                      NumericVector rowmeansA, NumericVector rowmeansB,
                      double meanA, double meanB) {
  // returns the dcov V-statistic kernel matrix H
  // for large bivariate samples (x, y)
  int i, j, k, m = I.size();
  NumericMatrix H(m, m);
  NumericVector xI(m), yI(m), aI(m), bI(m);
  
  for (i=0; i<m; i++) {
    k = I[i] - 1;
    xI(i) = x[k];
    yI(i) = y[k];
    aI(i) = rowmeansA[k];
    bI(i) = rowmeansB[k];
  }
  
  for (i=0; i<m; i++) {
    for (j=i; j<m; j++) {
      H(i, j) = ((fabs(xI(i)-xI(j)) - aI(i) - aI(j) + meanA) *
        (fabs(yI(i) - yI(j)) - bI(i) - bI(j) + meanB)) / ((double) m);
      H(j, i) = H(i, j);
    }
  }
  
  return H; 
}
