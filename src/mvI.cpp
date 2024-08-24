#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix Dxi(NumericMatrix Dx, IntegerVector ix);
  
// [[Rcpp::export]]
double Istat(NumericMatrix Dx, NumericMatrix Dy) {
  // compute independence coefficient I_n (the square root)
  // Dx and Dy are the Euclidean distance matrices
  
  int n = Dx.nrow();
  int i, j, k, m;
  double n2 = n*n, n3 = n*n2, n4 = n2*n2;
  double Cx, Cy, zd, zbar, z;
  IntegerVector ix(n), iy(n);
  NumericMatrix Dx2(n, n), Dy2(n, n);
  
  Cx = 0.0; Cy = 0.0; z = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      Cx += 2*Dx(i, j);
      Cy += 2*Dy(i, j);
      Dx2(i, j) = Dx(i, j) * Dx(i, j);
      Dy2(i, j) = Dy(i, j) * Dy(i, j);
      Dx2(j, i) = Dx2(i, j);
      Dy2(j, i) = Dy2(i, j);
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
  return sqrt((2 * zbar - zd - z) / (Cx + Cy - z));
}

NumericMatrix Dxi(NumericMatrix Dx, IntegerVector ix) {
  int i, j, n=Dx.nrow();
  NumericMatrix D(n, n);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      D(i, j) = Dx(ix[i], ix[j]);      
    }
  }
  return D;
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
  IntegerVector ix(n);
  NumericVector istats(R+1);
  NumericMatrix Dx2(n, n), Dy2(n, n), Dx2b(n, n);
  
  Cx = 0.0; Cy = 0.0; z = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      Cx += 2*Dx(i, j);
      Cy += 2*Dy(i, j);
      Dx2(i, j) = Dx(i, j) * Dx(i, j);
      Dy2(i, j) = Dy(i, j) * Dy(i, j);
      Dx2(j, i) = Dx2(i, j);
      Dy2(j, i) = Dy2(i, j);
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
    ix = sample(n, n, false, R_NilValue, false);
    Dx2b = Dxi(Dx2, ix);
    zd = 0.0; zbar = 0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        zd += sqrt(Dx2b(i, j) + Dy2(i, j));
        for (k=0; k<n; k++) {
          zbar += sqrt(Dx2b(k, i) + Dy2(k, j));
        }
      }
    }
    zd /= n2;
    zbar /= n3;
    istats[b] = sqrt((2 * zbar - zd - z) / (Cx + Cy - z));
  }
  return istats;
}

