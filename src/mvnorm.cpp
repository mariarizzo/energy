// Maria L. Rizzo <mrizzo@bgsu.edu>

#include <Rcpp.h>
using namespace Rcpp;

double   mvnEstat(NumericMatrix y);
double   sumdist(NumericMatrix x);

// M_ are constants defined in gnu c math.h

//[[Rcpp::export]]
double mvnEstat(NumericMatrix y)  {
  // compute E test statistic for multivariate normality
  // y is a *standardized* multivariate sample
  int    d = y.ncol(), n = y.nrow();
  int    i, j, k, maxterms=2000;
  double meanyy, meanyz, meanzz, stat;
  double delta, eps = 1.0e-7;
  double normy, yy, dif, sum, sum0, term;
  double D = (double) d;
  double lg0 = R::lgammafn(D / (double) 2);
  double lg1 = R::lgammafn((D+1.0) / (double) 2);
  double kd, logak, loggk;

  meanzz = 2.0 * exp(lg1 - lg0);  // the second mean

  // computing the first mean as series
  meanyz = 0.0;
  for (i=0; i<n; i++) {
      yy = 0.0;
      for (j=0; j<d; j++) {
        dif = y(i, j) * y(i, j);
        yy += dif;
      }
      normy = sqrt(yy);
      delta = 1.0;
      sum = 0.0;
      k = 0;
      while (delta > eps && k < maxterms) {
        kd = (double) k;
        sum0 = sum;
        logak = (kd+1) * log(yy) - R::lgammafn(kd+1) - kd*M_LN2 -
          log(2*kd+1) - log(2*kd+2);
        loggk = lg1 + R::lgammafn(kd+1.5) - R::lgammafn(kd+D/2+1);
        term = exp(logak + loggk);
        if (k % 2 == 0)
          sum += term;
          else
            sum -= term;
            delta = fabs(sum - sum0);
            k++;
      }
  if (delta < eps)
        meanyz += meanzz/M_SQRT2 + M_SQRT_2dPI * sum;
    else {
        meanyz += normy;
        Rf_warning("E|y-Z| did not converge, replaced by %f", normy);
    }
  }
  meanyz /= (double) n;

  meanyy = sumdist(y);  // computing third mean
  meanyy = (2.0 * meanyy / (double)(n*n));

  stat = ((double) n)*(2.0 * meanyz - meanzz - meanyy);
  return stat;
}


double sumdist(NumericMatrix x)
{
  // sum the pairwise distances between rows of data matrix x
  // without storing the distance matrix
  // lower triangle only
  // result is equivalent to this in R:  sum(dist(x))

  int n = x.nrow(), d = x.ncol();
  double s = 0.0, dsum, dif;
  for (int i=1; i<n; i++) {
    for (int j=0; j<i; j++) {
      dsum = 0.0;
      for (int k=0; k<d; k++) {
        dif = x(i, k) - x(j, k);
        dsum += dif * dif;
      }
      s += sqrt(dsum);
    }
  }
  return s;
}
