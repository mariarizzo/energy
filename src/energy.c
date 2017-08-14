/*
   energy.c: energy package

   Author:  Maria Rizzo <mrizzo @ bgnet.bgsu.edu>
   Created: 4 Jan 2004
   Updated: 2 April 2008    some functions moved to utilities.c
   Updated: 25 August 2016  mvnEstat converted to c++ in mvnorm.cpp

   ksampleEtest() performs the multivariate E-test for equal distributions,
                  complete version, from data matrix
   E2sample()     computes the 2-sample E-statistic without creating distance
   poisMstat()    computes the mean distance test of Poissonity
*/

#include <R.h>
#include <Rmath.h>

void   poisMstat(int *x, int *nx, double *stat);
void   ksampleEtest(double *x, int *byrow, int *nsamples, int *sizes, int *dim,
            int *R, double *e0, double *e, double *pval);
void   E2sample(double *x, int *sizes, int *dim, double *stat);

double edist(double **D, int m, int n);
double multisampleE(double **D, int nsamples, int *sizes, int *perm);
double twosampleE(double **D, int m, int n, int *xrows, int *yrows);
double E2(double **x, int *sizes, int *start, int ncol, int *perm);
double Eksample(double *x, int *byrow, int r, int d, int K, int *sizes, int *ix);
void   distance(double **bxy, double **D, int N, int d);

/* utilities.c */
extern double **alloc_matrix(int r, int c);
extern int    **alloc_int_matrix(int r, int c);
extern void   free_matrix(double **matrix, int r, int c);
extern void   free_int_matrix(int **matrix, int r, int c);

extern void   permute(int *J, int n);
extern void   roworder(double *x, int *byrow, int r, int c);
extern void   vector2matrix(double *x, double **y, int N, int d, int isroworder);

extern void   distance(double **bxy, double **D, int N, int d);
extern void   Euclidean_distance(double *x, double **Dx, int n, int d);
extern void   index_distance(double *x, double **Dx, int n, int d, double index);
extern void   sumdist(double *x, int *byrow, int *nrow, int *ncol, double *lowersum);


void poisMstat(int *x, int *nx, double *stat)
{
    /* computes the Poisson mean distance statistic */
    int i, j, k, n=(*nx);
    double eps=1.0e-10;
    double cvm, d, lambda, m, q;
    double Mcdf1, Mcdf0, Mpdf1, cdf1, cdf0;

    lambda = 0;
    for (i=0; i<n; i++)
        lambda += x[i];
    lambda /= ((double) n);
    q = qpois(1.0-eps, lambda, TRUE, FALSE) + 1;

    m = 0.0;
    for (j=0; j<n; j++) m += abs(x[j] - 1);
    m /= ((double) n);                   /* est of m_1 = E|1 - X| */
    Mcdf0 = (m + 1.0 - lambda) / 2.0;    /* M-est of F(0) */

    cdf0 = exp(-lambda);                 /* MLE of F(0) */
    d = Mcdf0 - cdf0;
    cvm = d * d * cdf0;   /* von Mises type of distance */

    for (i=1; i<q; i++) {
        m = 0;
        k = i + 1;
        for (j=0; j<n; j++) m += abs(x[j]-k);
        m /= ((double) n);  /* est of m_{i+1} = E|i+1 - X| */

        /* compute M-estimate of f(i) and F(i) */
        Mpdf1 = (m-(k-lambda)*(2.0*Mcdf0-1.0))/((double) 2.0*k);
        if (Mpdf1 < 0.0) Mpdf1 = 0.0;
        Mcdf1 = Mcdf0 + Mpdf1;
        if (Mcdf1 > 1) Mcdf1 = 1.0;

        cdf1 = ppois(i, lambda, TRUE, FALSE); /* MLE of F(i) */
        d = Mcdf1 - cdf1;
        cvm += d * d * (cdf1 - cdf0);

        cdf0 = cdf1;
        Mcdf0 = Mcdf1;
    }
    cvm *= n;
    *stat = cvm;
}


void E2sample(double *x, int *sizes, int *dim, double *stat) {
    /*
      compute test statistic *stat for testing H:F=G
      does not store distance matrix
      x must be in row order: x=as.double(t(x)) where
      x is pooled sample in matrix sum(en) by dim
    */
    int    m=sizes[0], n=sizes[1], d=(*dim);
    int    i, j, k, p, q;
    double dif, dsum, sumxx, sumxy, sumyy, w;

    sumxy = 0.0;
    for (i=0; i<m; i++) {
        p = i*d;
        for (j=m; j<m+n; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            sumxy += sqrt(dsum);
        }
    }
    sumxy /= (double)(m*n);
    sumxx = 0.0;
    for (i=1; i<m; i++) {
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            sumxx += sqrt(dsum);
        }
    }
    sumxx /= (double)(m*m);  /* half the sum */
    sumyy = 0.0;
    for (i=m+1; i<m+n; i++) {
        p = i*d;
        for (j=m; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            sumyy += sqrt(dsum);
        }
    }
    sumyy /= (double)(n*n);  /* half the sum */
    w = (double)(m*n)/(double)(m+n);
    *stat = 2.0*w*(sumxy - sumxx - sumyy);
}

void ksampleEtest(double *x, int *byrow,
                  int *nsamples, int *sizes, int *dim,
                  int *R, double *e0, double *e, double *pval)
{
    /*
      exported for R energy package: E test for equal distributions
      x         the pooled sample (or distances)
      *byrow    logical, TRUE if x is stored by row
                pass x=as.double(t(x)) if *byrow==TRUE
      *nsamples number of samples
      *sizes    vector of sample sizes
      *dim      dimension of data in x (0 if x is distance matrix)
      *R        number of replicates for permutation test
      *e0       observed E test statistic
      e         vector of replicates of E statistic
      *pval     approximate p-value
    */

    int    b, ek, i, k;
    int    B = (*R), K = (*nsamples), d=(*dim), N;
    int    *perm;
    double **data, **D;

    N = 0;
    for (k=0; k<K; k++)
        N += sizes[k];
    perm = Calloc(N, int);
    for (i=0; i<N; i++)
        perm[i] = i;
    D   = alloc_matrix(N, N);      /* distance matrix */
    if (d > 0) {
        data = alloc_matrix(N, d); /* sample matrix */
        vector2matrix(x, data, N, d, *byrow);
        distance(data, D, N, d);
        free_matrix(data, N, d);
    }
    else
        vector2matrix(x, D, N, N, *byrow);

    *e0 = multisampleE(D, K, sizes, perm);

    /* bootstrap */
    if (B > 0) {
        ek = 0;
        GetRNGstate();
        for (b=0; b<B; b++) {
            permute(perm, N);
            e[b] = multisampleE(D, K, sizes, perm);
            if ((*e0) < e[b]) ek++;
        }
        PutRNGstate();
        (*pval) = ((double) (ek + 1)) / ((double) (B + 1));
    }

    free_matrix(D, N, N);
    Free(perm);
}



double E2(double **x, int *sizes, int *start, int ncol, int *perm)
{
    int    m=sizes[0], n=sizes[1];
    int    row1=start[0], row2=start[1];
    int    i, j, k, p, q;
    double dif, dsum, sumxx, sumxy, sumyy, w;

    sumxy = 0.0;
    for (i=0; i<m; i++) {
        p = perm[row1 + i];
        for (j=0; j<n; j++) {
            dsum = 0.0;
            q = perm[row2 + j];
            for (k=0; k<ncol; k++) {
                dif = x[p][k] - x[q][k];
                dsum += dif * dif;
            }
            sumxy += sqrt(dsum);
        }
    }
    sumxy /= (double)(m * n);
    sumxx = 0.0;
    for (i=1; i<m; i++) {
        p = perm[row1 + i];
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = perm[row1 + j];
            for (k=0; k<ncol; k++) {
                dif = x[p][k] - x[q][k];
                dsum += dif * dif;
            }
            sumxx += sqrt(dsum);
        }
    }
    sumxx /= (double)(m * m);  /* half the sum */
    sumyy = 0.0;
    for (i=1; i<n; i++) {
        p = perm[row2 + i];
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = perm[row2 + j];
            for (k=0; k<ncol; k++) {
                dif = x[p][k] - x[q][k];
                dsum += dif * dif;
            }
            sumyy += sqrt(dsum);
        }
    }
    sumyy /= (double)(n * n);  /* half the sum */
    w = (double)(m * n)/(double)(m + n);
    return 2.0 * w * (sumxy - sumxx - sumyy);
}


double multisampleE(double **D, int nsamples, int *sizes, int *perm)
{
    /*
      returns the multisample E statistic
      D is square Euclidean distance matrix
      perm is a permutation of the row indices
    */
    int i, j, k, m, n;
    int *M;
    double e;

    M = Calloc(nsamples, int);
    M[0] = 0;
    for (k=1; k<nsamples; k++)
        M[k] = M[k-1] + sizes[k-1]; /* index where sample k begins */

    e = 0.0;
    for (i=0; i<nsamples; i++) {
        m = sizes[i];
        for (j=i+1; j<nsamples; j++) {
            n = sizes[j];
            e += twosampleE(D, m, n, perm+M[i], perm+M[j]);
        }
    }
    Free(M);
    return(e);
}

double twosampleE(double **D, int m, int n, int *xrows, int *yrows)
{
    /*
       return the e-distance between two samples
       corresponding to samples indexed xrows[] and yrows[]
       D is square Euclidean distance matrix
    */
    int    i, j;
    double sumxx=0.0, sumyy=0.0, sumxy=0.0;

    if (m < 1 || n < 1) return 0.0;
    for (i=0; i<m; i++)
        for (j=i+1; j<m; j++)
            sumxx += D[xrows[i]][xrows[j]];
    sumxx *= 2.0/((double)(m*m));
    for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
            sumyy += D[yrows[i]][yrows[j]];
    sumyy *= 2.0/((double)(n*n));
    for (i=0; i<m; i++)
        for (j=0; j<n; j++)
            sumxy += D[xrows[i]][yrows[j]];
    sumxy /= ((double) (m*n));

    return (double)(m*n)/((double)(m+n)) * (2*sumxy - sumxx - sumyy);
}

double edist(double **D, int m, int n)
{
    /*
      return the e-distance between two samples size m and n
      D is square Euclidean distance matrix
    */
    int    i, j;
    double sumxx=0.0, sumyy=0.0, sumxy=0.0;

    if (m < 1 || n < 1) return 0.0;
    for (i=0; i<m; i++)
        for (j=i+1; j<m; j++)
            sumxx += D[i][j];
    sumxx *= 2.0/((double)(m*m));
    for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
            sumyy += D[i][j];
    sumyy *= 2.0/((double)(n*n));
    for (i=0; i<m; i++)
        for (j=0; j<n; j++)
            sumxy += D[i][j];
    sumxy /= ((double) (m*n));
    return (double)(m*n)/((double)(m+n)) * (2*sumxy - sumxx - sumyy);
}

