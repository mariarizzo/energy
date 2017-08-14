/*
   utilities.c: some utilities for the energy package

   Author: Maria L. Rizzo
   (see energy package on CRAN
   or at personal.bgsu.edu/~mrizzo)

   alloc_matrix, alloc_int_matrix, free_matrix, free_int_matrix:
     use R (Calloc, Free) instead of C (calloc, free) for memory management

   permute            permutes the first n elements of an integer vector
   row_order          converts arg from column order to row order
   vector2matrix      copies double* arg into double** arg
   distance           computes Euclidean distance matrix from double**
   Euclidean_distance computes Euclidean distance matrix from double*
   index_distance     computes Euclidean distance matrix D then D^index
   sumdist            sums the distance matrix without creating the matrix

   Notes:
   1. index_distance (declaration and body of the function) revised in
      energy 1.3-0, 2/2011.
*/

#include <R.h>
#include <Rmath.h>

double **alloc_matrix(int r, int c);
int    **alloc_int_matrix(int r, int c);
void   free_matrix(double **matrix, int r, int c);
void   free_int_matrix(int **matrix, int r, int c);

void   permute(int *J, int n);
void   permute_check(int *J, int *N);
void   roworder(double *x, int *byrow, int r, int c);
void   vector2matrix(double *x, double **y, int N, int d, int isroworder);

void   distance(double **bxy, double **D, int N, int d);
void   Euclidean_distance(double *x, double **Dx, int n, int d);
void   index_distance(double **Dx, int n, double index);
void   sumdist(double *x, int *byrow, int *nrow, int *ncol, double *lowersum);



double **alloc_matrix(int r, int c)
{
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = Calloc(r, double *);
    for (i = 0; i < r; i++)
    matrix[i] = Calloc(c, double);
    return matrix;
}


int **alloc_int_matrix(int r, int c)
{
    /* allocate an integer matrix with r rows and c columns */
    int i;
    int **matrix;
    matrix = Calloc(r, int *);
    for (i = 0; i < r; i++)
    matrix[i] = Calloc(c, int);
    return matrix;
}

void free_matrix(double **matrix, int r, int c)
{
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) Free(matrix[i]);
    Free(matrix);
}

void free_int_matrix(int **matrix, int r, int c)
{
    /* free an integer matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) Free(matrix[i]);
    Free(matrix);
}

void permute(int *J, int n)
{
    /*
       permute the first n integers of J
       if n is length(J), equivalent to R:
            J <- rev(sample(J, length(J), replace=FALSE))
    */
    int i, j, j0, m=n;
    for (i=0; i<n-1; i++) {
		j = floor(runif(0.0, (double) m));
        m--;
        j0 = J[j];
        J[j] = J[m];
        J[m] = j0;
    }
}

void permute_check(int *J, int *N)
{
    /*
	the permute function can be called from R using this
	wrapper - the purpose is to check on a reported bug
    */
	int n = *N;
	permute(J, n);
}

void vector2matrix(double *x, double **y, int N, int d, int isroworder) {
    /* copy a d-variate sample into a matrix, N samples in rows */
    int i, k;
    if (isroworder == TRUE) {
        for (k=0; k<d; k++)
            for (i=0; i<N; i++)
                y[i][k] = (*(x+i*d+k));
        }
    else {
        for (k=0; k<N; k++)
            for (i=0; i<d; i++)
                y[i][k] = (*(x+k*N+i));
        }
    return;
}


void roworder(double *x, int *byrow, int r, int c) {
    /*
      utility to convert a vector from column order to row order
      assume that x is r by c matrix as a vector in column order
    */
    int    i, j, k, n=r*c;
    double *y;
    if (*byrow == TRUE) return;
    y = Calloc(n, double);
    i = 0;
    for (j=0; j<r; j++) {
        for (k=0; k<n; k+=r) {
            y[i] = x[k+j];
            i++;
        }
    }
    for (i=0; i<n; i++)
        x[i] = y[i];
    Free(y);
    *byrow = TRUE;
    return;
}


void distance(double **data, double **D, int N, int d) {
    /*
       compute the distance matrix of sample in N by d matrix data
       equivalent R code is:  D <- as.matrix(dist(data))
    */
    int    i, j, k;
    double dif;
    for (i=0; i<N; i++) {
        D[i][i] = 0.0;
        for (j=i+1; j<N; j++) {
            D[i][j] = 0.0;
            for (k=0; k<d; k++) {
                dif = data[i][k] - data[j][k];
                D[i][j] += dif*dif;
            }
            D[i][j] = sqrt(D[i][j]);
            D[j][i] = D[i][j];
        }
    }
    return;
}


void Euclidean_distance(double *x, double **Dx, int n, int d)
{
    /*
        interpret x as an n by d matrix, in row order (n vectors in R^d)
        compute the Euclidean distance matrix Dx
    */
    int i, j, k, p, q;
    double dsum, dif;
    for (i=1; i<n; i++) {
        Dx[i][i] = 0.0;
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            Dx[i][j] = Dx[j][i] = sqrt(dsum);
        }
    }
}


void index_distance(double **Dx, int n, double index)
{
    /*
        Dx is an n by n Euclidean distance matrix
        if index NEQ 1, compute D^index
    */
    int i, j;

    if (fabs(index - 1) > DBL_EPSILON) {
        for (i=0; i<n; i++)
            for (j=i+1; j<n; j++) {
                Dx[i][j] = R_pow(Dx[i][j], index);
                Dx[j][i] = Dx[i][j];
            }
    }
}


void sumdist(double *x, int *byrow, int *nrow, int *ncol, double *lowersum)
{
    /*
       sum all pairwise distances between rows of x
       equivalent to this in R:  h <- sum(dist(x))
       x must be in row order: x=as.double(t(x))
    */

    int i, j, k, p, q, n=(*nrow), d=(*ncol);
    double sum, dsum, dif;
    if (*byrow == FALSE)
        roworder(x, byrow, n, d);
    sum = 0.0;
    for (i=1; i<n; i++) {
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            sum += sqrt(dsum);
        }
    }
    (*lowersum) = sum;
}

