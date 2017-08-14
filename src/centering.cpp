// double centering utilities for the energy package
//
// Maria L. Rizzo <mrizzo@bgsu.edu>
// August, 2016



#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix D_center(NumericMatrix Dx);
NumericMatrix U_center(NumericMatrix Dx);

// [[Rcpp::export]]
NumericMatrix D_center(NumericMatrix Dx) {
  /*
  computes the double centered distance matrix for distance matrix Dx
   for dCov, dCor, etc.
  a_{ij} - a_{i.}/n - a_{.j}/n + a_{..}/n^2, all i, j
  */
  int j, k;
  int n = Dx.nrow();
  NumericVector akbar(n);
  NumericMatrix A(n, n);
  double abar = 0.0;

  for (k=0; k<n; k++) {
    akbar(k) = 0.0;
    for (j=0; j<n; j++) {
      akbar(k) += Dx(k, j);
    }
    abar += akbar(k);
    akbar(k) /= (double) (n);
  }
  abar /= (double) (n * n);

  for (k=0; k<n; k++)
    for (j=k; j<n; j++) {
      A(k, j) = Dx(k, j) - akbar(k) - akbar(j) + abar;
      A(j, k) = A(k, j);
    }

  return A;
}

// [[Rcpp::export]]
NumericMatrix U_center(NumericMatrix Dx) {
  /*
  computes the A_{kl}^U distances from the distance matrix (Dx_{kl}) for dCov^U
  U-centering: if Dx = (a_{ij}) then compute U-centered A^U using
  a_{ij} - a_{i.}/(n-2) - a_{.j}/(n-2) + a_{..}/((n-1)(n-2)), i \neq j
  and zero diagonal
  */
  int j, k;
  int n = Dx.nrow();
  NumericVector akbar(n);
  NumericMatrix A(n, n);
  double abar = 0.0;

  for (k=0; k<n; k++) {
    akbar(k) = 0.0;
    for (j=0; j<n; j++) {
      akbar(k) += Dx(k, j);
    }
    abar += akbar(k);
    akbar(k) /= (double) (n-2);
  }
  abar /= (double) ((n-1)*(n-2));

  for (k=0; k<n; k++)
    for (j=k; j<n; j++) {
      A(k, j) = Dx(k, j) - akbar(k) - akbar(j) + abar;
      A(j, k) = A(k, j);
    }
    /* diagonal is zero */
    for (k=0; k<n; k++)
      A(k, k) = 0.0;

  return A;
}

