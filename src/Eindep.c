/*
   Eindep.c: energy package
   E-statistics and test for multivariate independence (coef. I_n)
   Distance covariance and distance correlation: see dcov.c

   Author:   Maria Rizzo <mrizzo @ bgnet.bgsu.edu>
   Created:  June 15, 2004  (development)
   Last Modified:  April 5, 2008
*/

#include <R.h>
#include <Rmath.h>

void   indepE(double *x, double *y, int *byrow, int *dims, double *Istat);
void   indepEtest(double *x, double *y, int *byrow, int *dims,
                double *Istat, double *reps, double *pval);

void   squared_distance(double *x, double **D, int n, int d);

extern double **alloc_matrix(int r, int c);
extern int    **alloc_int_matrix(int r, int c);
extern void   free_matrix(double **matrix, int r, int c);
extern void   free_int_matrix(int **matrix, int r, int c);
extern void   permute(int *J, int n);
extern void   roworder(double *x, int *byrow, int r, int c);
extern void   Euclidean_distance(double *x, double **D, int n, int d);

void indepE(double *x, double *y, int *byrow, int *dims, double *Istat)
{
    /*
        E statistic for multiv. indep. of X in R^p and Y in R^q
        statistic returned is I_n^2 [nI_n^2 has a limit dist under indep]
        dims[0] = n (sample size)
        dims[1] = p (dimension of X)
        dims[2] = q (dimension of Y)
        Istat : the statistic I_n (normalized)
     */
    int    i, j, k, m, n, p, q;
    double Cx, Cy, Cz, C3, C4, n2, n3, n4, v;
    double **D2x, **D2y;

    n = dims[0];
    p = dims[1];
    q = dims[2];

    if (*byrow == FALSE) {
        /* avoid this step: use as.double(t(x)) in R */
        roworder(x, byrow, n, p);
        *byrow = FALSE;  /* false for y */
        roworder(y, byrow, n, q);
    }

    D2x = alloc_matrix(n, n);
    D2y = alloc_matrix(n, n);
    Euclidean_distance(x, D2x, n, p);
    Euclidean_distance(y, D2y, n, q);

    Cx = Cy = Cz = C3 = C4 = 0.0;
    n2 = ((double) n) * n;
    n3 = n2 * n;
    n4 = n2 * n2;

    /* compute observed test statistic */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            Cx += (D2x[i][j]);
            Cy += (D2y[i][j]);
            Cz += sqrt(D2x[i][j]*D2x[i][j] + D2y[i][j]*D2y[i][j]);
        }
    }
    Cx = 2.0 * Cx / n2;
    Cy = 2.0 * Cy / n2;
    Cz = 2.0 * Cz / n2;

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++) {
                C3 += sqrt(D2x[k][i]*D2x[k][i] + D2y[k][j]*D2y[k][j]);
                for (m=0; m<n; m++)
                    C4 += sqrt(D2x[i][k]*D2x[i][k] + D2y[j][m]*D2y[j][m]);
            }
        }
    }
    C3 /= n3;
    C4 /= n4;
    v = Cx + Cy - C4;
    *Istat = (2.0 * C3 - Cz - C4) / v;
    free_matrix(D2x, n, n);
    free_matrix(D2y, n, n);
    return;
}


void indepEtest(double *x, double *y, int *byrow, int *dims,
                double *Istat, double *reps, double *pval) {
    /*
        approx permutation E test for multiv. indep. of X in R^p and Y in R^q
        statistic is I_n^2, where nI_n^2 -> Q
        dims[0] = n (sample size)
        dims[1] = p (dimension of X)
        dims[2] = q (dimension of Y)
        dims[3] = B (number of replicates, dimension of reps)
        Istat : the statistic I_n (normalized)
     */
    int    b, i, j, k, m, n, p, q, B, M;
    int    *perm;
    double Cx, Cy, Cz, C3, C4, n2, n3, n4, v;
    double **D2x, **D2y;

    n = dims[0];
    p = dims[1];
    q = dims[2];
    B = dims[3];

    if (*byrow == FALSE) {
        /* avoid this step: use as.double(t(x)) in R */
        roworder(x, byrow, n, p);
        *byrow = FALSE;  /* false for y */
        roworder(y, byrow, n, q);
    }

    D2x = alloc_matrix(n, n);
    D2y = alloc_matrix(n, n);
    squared_distance(x, D2x, n, p);
    squared_distance(y, D2y, n, q);

    Cx = Cy = Cz = C3 = C4 = 0.0;
    n2 = ((double) n) * n;
    n3 = n2 * n;
    n4 = n2 * n2;

    /* compute observed test statistic */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            Cx += sqrt(D2x[i][j]);
            Cy += sqrt(D2y[i][j]);
            Cz += sqrt(D2x[i][j] + D2y[i][j]);
        }
    }
    Cx = 2.0 * Cx / n2;
    Cy = 2.0 * Cy / n2;
    Cz = 2.0 * Cz / n2;

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++) {
                C3 += sqrt(D2x[k][i] + D2y[k][j]);
                for (m=0; m<n; m++)
                    C4 += sqrt(D2x[i][k] + D2y[j][m]);
            }
        }
    }
    C3 /= n3;
    C4 /= n4;

    v = Cx + Cy - C4;
    *Istat = (2.0 * C3 - Cz - C4) / v;

    M = 0;
    /* compute the replicates */
    if (B > 0) {
        GetRNGstate();
        perm = Calloc(n, int);
        for (i=0; i<n; i++)
            perm[i] = i;
        for (b = 0; b < B; b++) {
            permute(perm, n);
            Cz = 0.0;
            C3 = 0.0;
            for (i=0; i<n; i++)
                for (j=0; j<n; j++) {
                    Cz += sqrt(D2x[i][j] + D2y[perm[i]][perm[j]]);
                    for (k=0; k<n; k++) {
                        C3 += sqrt(D2x[k][perm[i]] + D2y[k][perm[j]]);
                    }
                }
            Cz /= n2;
            C3 /= n3;
            reps[b] = (2.0 * C3 - Cz - C4) / v;
            if (reps[b] >= (*Istat)) M++;
        }
        *pval = (double) M / (double) B;
        PutRNGstate();
        Free(perm);
    }

    free_matrix(D2x, n, n);
    free_matrix(D2y, n, n);
    return;
}


void squared_distance(double *x, double **D2, int n, int d)
{
    /*
        interpret x as an n by d matrix, in row order (n vectors in R^d)
        compute the squared distance matrix D2
    */
    int i, j, k, p, q;
    double dsum, dif;
    for (i=1; i<n; i++) {
        D2[i][i] = 0.0;
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            D2[i][j] = D2[j][i] = dsum;
        }
    }
}

