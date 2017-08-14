#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix U_center(NumericMatrix);
double        U_product(NumericMatrix, NumericMatrix);

// [[Rcpp::export]]
NumericMatrix projection(NumericMatrix Dx, NumericMatrix Dz) {
  /*
  returns the projection of A(x) distance matrix Dx onto the
  orthogonal complement of C(z) distance matrix;
  both Dx and Dz are n by n distance or dissimilarity matrices
  the projection is an n by n matrix
  */
  int    n = Dx.nrow();
  int    i, j;
  NumericMatrix A(n, n), C(n, n), P(n, n);
  double AC, CC, c1;
  double eps = std::numeric_limits<double>::epsilon();  //machine epsilon

  A = U_center(Dx);        // U-centering to get A^U etc.
  C = U_center(Dz);
  AC = U_product(A, C);    // (A,C) = dcov^U
  CC = U_product(C, C);
  c1 = 0.0;
  // if (C,C)==0 then C==0 so c1=(A,C)=0
  if (fabs(CC) > eps)
    c1 = AC / CC;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      P(i, j) = A(i, j) - c1 * C(i, j);
    }
    return P;
}

