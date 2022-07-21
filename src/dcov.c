/*
   dcov.c: distance correlation and covariance statistics
   and dCov test for multivariate independence

   Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007)
   "Measuring and testing dependence by correlation of distances"
   Annals of Statistics, Vol. 35 No. 6, pp. 2769-2794.

 Author: Maria L. Rizzo
 energy package
 github.com/mariarizzo/energy
 
 */


#include <R.h>
#include <Rmath.h>

void   dCOVtest(double *x, double *y, int *nrow, int *nreps,
                double *reps, double *DCOV, double *pval);
void   dCOV(double *x, double *y, int *nrow, double *DCOV);
double Akl(double **akl, double **A, int n);

/* functions in utilities.c */
extern double **alloc_matrix(int r, int c);
extern int    **alloc_int_matrix(int r, int c);
extern void   free_matrix(double **matrix, int r, int c);
extern void   free_int_matrix(int **matrix, int r, int c);
extern void   permute(int *J, int n);
extern void   roworder(double *x, int *byrow, int r, int c);
extern void   Euclidean_distance(double *x, double **Dx, int n, int d);
extern void   index_distance(double **Dx, int n, double index);
extern void   vector2matrix(double *x, double **y, int N, int d, int isroworder);


void dCOVtest(double *x, double *y, int *nrow, int *nreps,
              double *reps, double *DCOV, double *pval) {
    /*  input vectors must expand to distance matrices 
        any exponent must be pre-computed in R
        computes dCov(x,y), dCor(x,y), dVar(x), dVar(y)
        V-statistic is n*dCov^2 where n*dCov^2 --> Q
        DCOV  : vector [dCov, dCor, dVar(x), dVar(y), mean(A), mean(B)]
     */
    int    i, j, k, r, J, K, M;
    int    n = nrow[0], R = nreps[0];
    int*   perm;
    double **Dx, **Dy, **A, **B;
    double dcov, V;
    double n2 = (double) n * n;
    Dx = alloc_matrix(n, n);
    Dy = alloc_matrix(n, n);
  	vector2matrix(x, Dx, n, n, 1);
	  vector2matrix(y, Dy, n, n, 1);

    A = alloc_matrix(n, n);
    B = alloc_matrix(n, n);
    Akl(Dx, A, n);
    Akl(Dy, B, n);
    free_matrix(Dx, n, n);
    free_matrix(Dy, n, n);

    /* compute dCov(x,y), dVar(x), dVar(y) */
    for (k=0; k<4; k++)
        DCOV[k] = 0.0;
    for (k=0; k<n; k++)
        for (j=0; j<n; j++) {
            DCOV[0] += A[k][j]*B[k][j];
            DCOV[2] += A[k][j]*A[k][j];
            DCOV[3] += B[k][j]*B[k][j];
        }

    for (k=0; k<4; k++) {
        DCOV[k] /= n2;
        if (DCOV[k] > 0)
            DCOV[k] = sqrt(DCOV[k]);
            else DCOV[k] = 0.0;
    }
    /* compute dCor(x, y) */
    V = DCOV[2]*DCOV[3];
    if (V > DBL_EPSILON)
        DCOV[1] = DCOV[0] / sqrt(V);
        else DCOV[1] = 0.0;

	if (R > 0) {
		/* compute the replicates */
		if (DCOV[1] > 0.0) {
			perm = Calloc(n, int);
			M = 0;
			for (i=0; i<n; i++) perm[i] = i;
      GetRNGstate();
      for (r=0; r<R; r++) {
			   permute(perm, n);
			   dcov = 0.0;
			   for (k=0; k<n; k++) {
				   K = perm[k];
				   for (j=0; j<n; j++) {
					   J = perm[j];
					   dcov += A[k][j]*B[K][J];
				   }
			   }
			   dcov /= n2;
			   dcov = sqrt(dcov);
			   reps[r] = dcov;
			   if (dcov >= DCOV[0]) M++;
			}
			*pval = (double) (M+1) / (double) (R+1);
      PutRNGstate();
			Free(perm);
		} else {
		    *pval = 1.0;
			}
	}

    free_matrix(A, n, n);
    free_matrix(B, n, n);
    return;
}

void dCOV(double *x, double *y, int *nrow, double *DCOV) {
  /*  input vectors must expand to distance matrices 
      any exponent must be pre-computed in R
        computes dCov(x,y), dCor(x,y), dVar(x), dVar(y)
        V-statistic is n*dCov^2 where n*dCov^2 --> Q
        DCOV  : vector [dCov, dCor, dVar(x), dVar(y)]
     */

    int    j, k, n = nrow[0];
    double **Dx, **Dy, **A, **B;
    double V, n2 = (double) n * n;

    Dx = alloc_matrix(n, n);
    Dy = alloc_matrix(n, n);
		vector2matrix(x, Dx, n, n, 1);
		vector2matrix(y, Dy, n, n, 1);

    A = alloc_matrix(n, n);
    B = alloc_matrix(n, n);
    Akl(Dx, A, n);
    Akl(Dy, B, n);
    free_matrix(Dx, n, n);
    free_matrix(Dy, n, n);

    n2 = ((double) n) * n;

    /* compute dCov(x,y), dVar(x), dVar(y) */
    for (k=0; k<4; k++)
        DCOV[k] = 0.0;
    for (k=0; k<n; k++)
        for (j=0; j<n; j++) {
            DCOV[0] += A[k][j]*B[k][j];
            DCOV[2] += A[k][j]*A[k][j];
            DCOV[3] += B[k][j]*B[k][j];
        }

    for (k=0; k<4; k++) {
        DCOV[k] /= n2;
        if (DCOV[k] > 0)
            DCOV[k] = sqrt(DCOV[k]);
            else DCOV[k] = 0.0;
    }
    /* compute dCor(x, y) */
    V = DCOV[2]*DCOV[3];
    if (V > DBL_EPSILON)
        DCOV[1] = DCOV[0] / sqrt(V);
        else DCOV[1] = 0.0;

    free_matrix(A, n, n);
    free_matrix(B, n, n);
    return;
}

double Akl(double **akl, double **A, int n) {
    /* -computes the A_{kl} or B_{kl} distances from the
        distance matrix (a_{kl}) or (b_{kl}) for dCov, dCor, dVar
        dCov = mean(Akl*Bkl), dVar(X) = mean(Akl^2), etc.
    */
    int j, k;
    double *akbar;
    double abar;

    akbar = Calloc(n, double);
    abar = 0.0;
    for (k=0; k<n; k++) {
        akbar[k] = 0.0;
        for (j=0; j<n; j++) {
            akbar[k] += akl[k][j];
        }
        abar += akbar[k];
        akbar[k] /= (double) n;
    }
    abar /= (double) (n*n);

    for (k=0; k<n; k++)
        for (j=k; j<n; j++) {
            A[k][j] = akl[k][j] - akbar[k] - akbar[j] + abar;
            A[j][k] = A[k][j];
        }
    Free(akbar);
    return(abar);
}
