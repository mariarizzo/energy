#include <Rcpp.h>
using namespace Rcpp;

int kgroups_update(NumericMatrix x, int k, IntegerVector clus,
                   IntegerVector sizes, NumericVector within, bool distance);
List kgroups_start(NumericMatrix x, int k, IntegerVector clus,
                   int iter_max, bool distance);

int kgroups_update(NumericMatrix x, int k, IntegerVector clus,
                      IntegerVector sizes, NumericVector w, bool distance) {
  /*
   * k-groups one pass through sample moving one point at a time
   * x: data matrix or distance
   * k: number of clusters
   * clus: clustering vector clus(i)==j ==> x_i is in cluster j
   * sizes: cluster sizes
   * within: vector of within cluster dispersions
   * distance: true if x is distance matrix
   * update clus, sizes, and withins
   * return count = number of points moved
   */

  int n = x.nrow(), d = x.ncol();
  int i, j, I, J, ix, nI, nJ;
  NumericVector rowdst(k), e(k);
  int best, count = 0;
  double dsum, dif;

  for (ix = 0; ix < n; ix++) {
    I = clus(ix);
    nI = sizes(I);
    if (nI > 1) {
      // calculate the E-distances of this point to each cluster
      rowdst.fill(0.0);
      for (i = 0; i < n; i++) {
        J = clus(i);
        if (distance == true) {
          rowdst(J) += x(ix, i);
        } else {
          dsum = 0.0;
          for (j = 0; j < d; j++) {
            dif = x(ix, j) - x(i, j);
            dsum += dif * dif;
          }
          rowdst(J) += sqrt(dsum);
        }
      }

      for (J = 0; J < k; J++) {
        nJ = sizes(J);
        e(J) = (2.0 / (double) nJ) * (rowdst(J) - w(J));
      }

      best = Rcpp::which_min(e);
      if (best != I) {
        // move this point and update
        nI = sizes(I);
        nJ = sizes(best);
        w(best) = (((double) nJ) * w(best) + rowdst(best)) / ((double) (nJ + 1));
        w(I) = (((double) nI) * w(I) - rowdst(I)) / ((double) (nI - 1));
        clus(ix) = best;
        sizes(I) = nI - 1;
        sizes(best) = nJ + 1;
        count ++;  // number of moves
        }
      }
    }

  return count;
}



// [[Rcpp::export]]
List kgroups_start(NumericMatrix x, int k, IntegerVector clus,
                   int iter_max, bool distance) {
  // k-groups clustering with initial clustering vector clus
  // up to iter_max iterations of n possible moves each
  // distance: true if x is distance matrix
    NumericVector within(k, 0.0);
  IntegerVector sizes(k, 0);
  double dif, dsum;
  int I, J, h, i, j;
  int n = x.nrow(), d = x.ncol();

  for (i = 0; i < n; i++) {
    I = clus(i);
    sizes(I)++;
    for (j = 0; j < i; j++) {
      J = clus(j);
      if (I == J) {
        if (distance == true) {
          within(I) += x(i, j);
        } else {
          dsum = 0.0;
          for (h = 0; h < d; h++) {
            dif = x(i, h) - x(j, h);
            dsum += dif * dif;
          }
          within(I) += sqrt(dsum);
        }
      }
    }
  }
  for (I = 0; I < k; I++)
    within(I) /= ((double) sizes(I));

  int it = 1, count = 1;
  count = kgroups_update(x, k, clus, sizes, within, distance);

  while (it < iter_max && count > 0) {
    count = kgroups_update(x, k, clus, sizes, within, distance);
    it++;
  }
  double W = Rcpp::sum(within);

  return List::create(
        _["within"] = within,
        _["W"] = W,
        _["sizes"] = sizes,
        _["cluster"] = clus,
        _["iterations"] = it,
        _["count"] = count);
}


