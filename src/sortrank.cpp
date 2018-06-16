#include <Rcpp.h>
using namespace Rcpp;

/*
 *  sort_rank() returns the order(), ranks() and sorted input vector
 *    in a Rcpp::List
 *  as a substitute for energy::sortrank() function in R
 *  benchmarks show the R implementation is faster
 *  order_rank() returns an IntegerMatrix of order(x), rank(x)
*/

List call_sortrank(NumericVector x, Function f);
IntegerMatrix order_rank(NumericVector x);
List sort_rank(NumericVector x);
    
List sort_rank(NumericVector x) {
  // will not print a warning about ties, may differ from R order()
  // almost as fast as my sortrank() function in R
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  IntegerVector o = match(y, x);
  IntegerVector N = seq(1, x.size());
  IntegerVector r = match(N, o);
  Rcpp::List L = Rcpp::List::create(
    _["x"] = y,
    _["ix"] = o,
    _["r"] = r);
  return L;
}

IntegerMatrix order_rank(NumericVector x) {
  // get the order(x) and rank(x) vectors
  // return in matrix [order, rank]
  int n = x.size();
  IntegerMatrix R(n, 2);
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  IntegerVector o = match(y, x);
  IntegerVector N = seq(1, x.size());
  IntegerVector r = match(N, o);
  R(_, 0) = o;
  R(_, 1) = r;
  return R;
}    

List call_sortrank(NumericVector x, Function f) {
  return f(x);
}


