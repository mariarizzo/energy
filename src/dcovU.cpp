#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix U_center(NumericMatrix);

//[[Rcpp::export]]
NumericVector dcovU_stats(NumericMatrix Dx, NumericMatrix Dy) {
  // x and y must be square distance matrices
  NumericMatrix A = U_center(Dx);
  NumericMatrix B = U_center(Dy);
  double ab = 0.0, aa = 0.0, bb = 0.0;
  double V, dcorU = 0.0;
  double eps = std::numeric_limits<double>::epsilon();  //machine epsilon
  int n = Dx.nrow();
  int n2 = n * (n - 3);

  for (int i=0; i<n; i++)
    for (int j=0; j<i; j++) {
      // U-centered is symmetric, with zero diagonal
      ab += A(i, j) * B(i, j);
      aa += A(i, j) * A(i, j);
      bb += B(i, j) * B(i, j);
    }
  ab = 2.0 * ab / (double) n2;
  aa = 2.0 * aa / (double) n2;
  bb = 2.0 * bb / (double) n2;
  V = aa * bb;
  if (V > eps)
    dcorU = ab / sqrt(V);

  return NumericVector::create(
    _["dCovU"] = ab,
    _["bcdcor"] = dcorU,
    _["dVarXU"] = aa,
    _["dVarYU"] = bb
  );
}

