#include <Rcpp.h>
using namespace Rcpp;

// alternate algorithms to compute dcov V-statistics or U-statistics
// without storing distance matrices 
// algebraically equivalent to dcor(x,y)^2, bcdcor, dcov(x,y)^2, dcovU(x,y), etc.


NumericVector Btree_sum (IntegerVector y, NumericVector z);

double sum_paired_dist(NumericMatrix x, NumericMatrix y);
double sum_dist3(NumericMatrix x, NumericMatrix y);
NumericVector rowSums_dist(NumericMatrix x);

NumericVector get_dcov_sums (NumericMatrix x, NumericMatrix y, bool all);
NumericVector dcov_UV(NumericMatrix x, NumericMatrix y, bool all);
  
// [[Rcpp::export]]
NumericVector dcov_UV(NumericMatrix x, NumericMatrix y, bool all=true) {
  // compute the dcov U statistic and the dcov V statistic
  // by alternate formula without storing distance matrices
  // returns the V-statistic(s) dcov_(x,y)^2
  // returns the U-statistic(s), unbiased for dcov^2(X,Y)
  // only the divisors differ for V and U formulas
  // optionally compute the distance variances, dcor^2 and bcdcor if all=true

  NumericVector stats =
    NumericVector::create(_["dcovV"]=0.0, _["dcorV"]=NA_REAL, 
                          _["Vx"]=NA_REAL, _["Vy"]=NA_REAL,
                          _["dcovU"]=0.0, _["dcorU"]=NA_REAL, 
                          _["Ux"]=NA_REAL, _["Uy"]=NA_REAL);
  int n = x.nrow();
  int d1, d2, d3, D1, D2, D3;
  d1 = n * n;
  d3 = d1 * n;
  d2 = d3 * n;
  D1 = n * (n - 3);
  D3 = D1 * (n - 2);
  D2 = D3 * (n - 1);
  
  NumericVector sums = get_dcov_sums(x, y, all);
  stats["dcovV"] = sums["sumAB"] / ((double) d1) + 
    sums["sumAsumB"] / ((double) d2) - 2.0 * sums["rowsAB"] / ((double) d3);
  stats["dcovU"] = sums["sumAB"] / ((double) D1) + 
    sums["sumAsumB"] / ((double) D2) - 2.0 * sums["rowsAB"] / ((double) D3);

  if (all == true) {
    // compute distance variances
    stats["Vx"] = sums["sumAA"] / ((double) d1) + 
      sums["sumAsumA"] / ((double) d2) - 2.0 * sums["rowsAA"] / ((double) d3);
    stats["Vy"] = sums["sumBB"] / ((double) d1) + 
      sums["sumBsumB"] / ((double) d2) - 2.0 * sums["rowsBB"] / ((double) d3);
    stats["Ux"] = sums["sumAA"] / ((double) D1) + 
      sums["sumAsumA"] / ((double) D2) - 2.0 * sums["rowsAA"] / ((double) D3);
    stats["Uy"] = sums["sumBB"] / ((double) D1) + 
      sums["sumBsumB"] / ((double) D2) - 2.0 * sums["rowsBB"] / ((double) D3);
    double dvars = stats["Vx"] * stats["Vy"];
    stats["dcorV"] = 0.0;
    if (dvars > 0) 
      stats["dcorV"] = stats["dcovV"] / sqrt(dvars);
    dvars = stats["Ux"] * stats["Uy"];
    stats["dcorU"] = 0.0;
    if (dvars > 0) 
      stats["dcorU"] = stats["dcovU"] / sqrt(dvars);
  }
  return stats;
}
  
// [[Rcpp::export]]
NumericVector get_dcov_sums (NumericMatrix x, NumericMatrix y,
                             bool all=true) {
  // get sums needed to compute dcov and dvar stats (V or U)
  // without storing distance matrices
  // V = sumAB/n^2 + sumAsumB/n^4 - 2/n^3 rowsAB
  // if all==FALSE, do not compute stats for distance variance
  NumericVector sums =
    NumericVector::create(_["sumAB"]=0.0, _["sumAsumB"]=0.0, _["rowsAB"]=0.0,
                          _["sumAA"]=NA_REAL, _["sumAsumA"]=NA_REAL, _["rowsAA"]=NA_REAL,
                          _["sumBB"]=NA_REAL, _["sumBsumB"]=NA_REAL, _["rowsBB"]=NA_REAL);
  double dcov_sum3 = 0.0, dvarx_sum3 = 0.0, dvary_sum3 = 0.0;
  int i, n = x.nrow();
  NumericVector rowsumsA = rowSums_dist(x);
  NumericVector rowsumsB = rowSums_dist(y);
  for (i = 0; i < n; i++) {
    dcov_sum3 += rowsumsA(i) * rowsumsB(i);
  }
  double sumA = sum(rowsumsA);
  double sumB = sum(rowsumsB);
  double sumAB = sum_paired_dist(x, y);

  if (all == true) {
    // optional computations if distance variance is needed
    sums["sumAA"] = sum_paired_dist(x, x);
    sums["sumBB"] = sum_paired_dist(y, y);
    for (int i = 0; i < n; i++) {
      dvarx_sum3 += rowsumsA(i) * rowsumsA(i);
      dvary_sum3 += rowsumsB(i) * rowsumsB(i);
    }
    sums["rowsAA"] = dvarx_sum3;
    sums["rowsBB"] = dvary_sum3;
    sums["sumAsumA"] = sumA * sumA;
    sums["sumBsumB"] = sumB * sumB;
  }

  sums["sumAB"] = sumAB;
  sums["sumAsumB"] = sumA * sumB;
  sums["rowsAB"] = dcov_sum3;
  return sums;
}

