#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calc_dist(NumericMatrix x) {
  int n = x.nrow(), d = x.ncol(), i, j, k;
  double dsum, dk;
  NumericMatrix Dx(n, n);
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      if (i == j) {
        Dx(i, i) = 0.0;
      } else {
        dsum = 0.0;
        for (k = 0; k < d; k++) {
          dk = x(i,k) - x(j,k);
          dsum += dk * dk;
        }
        Dx(i, j) = sqrt(dsum);
        Dx(j, i) = sqrt(dsum);
      }
    }
  }
  return Dx;
}

