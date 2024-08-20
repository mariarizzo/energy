#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Istat(NumericMatrix Dx, NumericMatrix Dy) {
  // compute independence coefficient I_n (the square root)
  // Dx and Dy are the Euclidean distance matrices
  int n = Dx.nrow();
  int i, j, k, m;
  double n2 = n*n, n3 = n*n2, n4 = n2*n2;
  double Cx, Cy, zd, zbar, z;
  NumericMatrix Dx2(n, n), Dy2(n, n);
  
  Cx = 0.0; Cy = 0.0; zd = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Cx += Dx(i, j);
      Cy += Dy(i, j);
      Dx2(i, j) = Dx(i, j) * Dx(i, j);
      Dy2(i, j) = Dy(i, j) * Dy(i, j);
      zd += sqrt(Dx2(i, j) + Dy2(i, j));
    }
  }
  Cx /= n2;
  Cy /= n2;
  zd /= n2;
  
  zbar = 0.0; z = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      for (k=0; k<n; k++) {
        zbar += sqrt(Dx2(k, i) + Dy2(k, j));
        for (m=0; m<n; m++) {
          z += sqrt(Dx2(k, i) + Dy2(m, j));
        }
      }
    }
  }
  zbar /= n3;
  z /= n4;
  return sqrt((2 * zbar - zd - z) / (Cx + Cy - z));
}


// [[Rcpp::export]]
NumericVector Istats(NumericMatrix Dx, NumericMatrix Dy, int R) {
// computes the observed I_n independence coefficient and
// R permutation replicates
// denominator is invariant to permutations of row indices
// returns vector c(I_n, I_n^*) 
  int n = Dx.nrow();
  int b, i, j, k, m;
  double n2 = n*n, n3 = n*n2, n4 = n2*n2;
  double Cx, Cy, zd, zbar, z;
  IntegerVector idx(n);
  NumericVector istats(R+1);
  NumericMatrix Dx2(n, n), Dy2(n, n);
  
  Cx = 0.0; Cy = 0.0; z = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Cx += Dx(i, j);
      Cy += Dy(i, j);
      Dx2(i, j) = Dx(i, j) * Dx(i, j);
      Dy2(i, j) = Dy(i, j) * Dy(i, j);
    }
  }
  Cx /= n2;
  Cy /= n2;
  
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      for (k=0; k<n; k++) {
        for (m=0; m<n; m++) {
          z += sqrt(Dx2(k, i) + Dy2(m, j));
        }
      }
    }
  }
  z /= n4;
  
  zbar = 0.0; zd = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      zd += sqrt(Dx2(i, j) + Dy2(i, j));
      for (k=0; k<n; k++) {
        zbar += sqrt(Dx2(k, i) + Dy2(k, j));
      }
    }
  }
  zbar /= n3;
  zd /= n2;
  istats[0] = sqrt((2 * zbar - zd - z) / (Cx + Cy - z));

  for (b=1; b<=R; b++) {
    int ii, jj, kk;
    idx = sample(n, n, false, R_NilValue, false);
    zbar = 0.0; zd = 0.0;
    for (i=0; i<n; i++) {
      ii = idx[i];
      for (j=0; j<n; j++) {
        jj = idx[j];
        zd += sqrt(Dx2(ii, jj) + Dy2(i, j));
        for (k=0; k<n; k++) {
          kk = idx[k];
          zbar += sqrt(Dx2(kk, ii) + Dy2(k, j));
        }
      }
    }
    zbar /= n3;
    zd /= n2;
    istats[b] = sqrt((2 * zbar - zd - z) / (Cx + Cy - z));
  }
  return istats;
}
