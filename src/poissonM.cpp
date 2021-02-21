#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(.poisMstat)]]
NumericVector poisMstat(IntegerVector x)
{
  /* computes the Poisson mean distance statistic */
  int i, j, k, n=x.size();
  double eps=1.0e-10;
  double ad, cvm, d, lambda, m, q;
  double Mcdf1, Mcdf0, Mpdf1, cdf1, cdf0;
  NumericVector stats(2);

  lambda = mean(x);
  q = R::qpois(1.0-eps, lambda, TRUE, FALSE) + 1;

  m = 0.0;
  for (j=0; j<n; j++) m += abs(x(j) - 1);
  m /= ((double) n);                   /* est of m_1 = E|1 - X| */
  Mcdf0 = (m + 1.0 - lambda) / 2.0;    /* M-est of F(0) */

  cdf0 = exp(-lambda);                 /* MLE of F(0) */
  d = Mcdf0 - cdf0;
  cvm = d * d * cdf0;   /* von Mises type of distance */
  ad = d * d * cdf0 / (cdf0 * (1-cdf0));  /* A-D weight */
    
  for (i=1; i<q; i++) {
    m = 0;
    k = i + 1;
    for (j=0; j<n; j++) m += abs(x(j)-k);
    m /= ((double) n);  /* est of m_{i+1} = E|i+1 - X| */

  /* compute M-estimate of f(i) and F(i) */
  Mpdf1 = (m-(k-lambda)*(2.0*Mcdf0-1.0))/((double) 2.0*k);
  if (Mpdf1 < 0.0) Mpdf1 = 0.0;
  Mcdf1 = Mcdf0 + Mpdf1;
  if (Mcdf1 > 1) Mcdf1 = 1.0;

  cdf1 = R::ppois(i, lambda, TRUE, FALSE); /* MLE of F(i) */
  d = Mcdf1 - cdf1;
  cvm += d * d * (cdf1 - cdf0);
  ad += d * d * (cdf1 - cdf0) / (cdf1 * (1-cdf1));
  
  cdf0 = cdf1;
  Mcdf0 = Mcdf1;
  }
  cvm *= n;
  ad *= n;
  stats(0) = cvm;
  stats(1) = ad;
  return stats;
}
