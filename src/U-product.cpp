#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double U_product(NumericMatrix U, NumericMatrix V) {
  // U and V are U-centered dissimilarity matrices of the two samples
  int n = U.nrow();
  int i, j;
  double sums = 0.0;

  for (i = 0; i < n; i++)
    for (j=0; j<i; j++)
      sums += U(i, j) * V(i, j);
  sums = 2.0 * sums / ((double) n * (n-3));
  return (sums);
}

