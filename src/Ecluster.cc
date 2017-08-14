/*
   Ecluster.cc: energy package

   Author:  Maria Rizzo <mrizzo @ bgnet.bgsu.edu>
   Created: Created: 12 Dec 2002
   Revised: 4 Jan 2004 for R-1.8.1 energy package
   Revised: 28 Jan 2004
*/

#include "R.h"
#include "Rmath.h"
#include "ECl.h"

extern "C" {

double **alloc_matrix(int r, int c);
void   free_matrix(double **matrix, int r, int c);
void   Emin_hclust(double *diss, int *en, int *merge, double *height, int *order);
void   lower2square(double **dst, double *diss, int n);

void Emin_hclust(double *diss, int *en, int *merge, double *height, int *order)
{
    // performs hierarchical E-clustering by minimum cluster E-distance
    //    diss    lower.tri of n by n distance matrix, column order
    //    en      sample size n
    //    merge   (n-1) by 2 array as a vector in col order; see hclust object
    //    height  vector length n-1; see hclust object
    //    order   vector length n; see hclust object

    int    i, step, I, J;
    int    n = (*en);
    double e;
    double *E;
    double **Edst;
    double **dst;
    int    *m1, *m2;
    ECl    c;  //clustering object

    c.init(n);

    dst = alloc_matrix(n, n);
    Edst = alloc_matrix(n, n); //E dist between clusters
    E =  Calloc(n, double);
    m1 = Calloc(n-1, int);
    m2 = Calloc(n-1, int);

    // convert lower.tri in vector form to square matrix
    lower2square(dst, diss, n);

    //E-hierarchical clustering

    E[0] = c.init_Edst(dst, Edst);
    step = 0;
    while (c.len() > 1)
    {
        e = c.merge_minEdst(dst, Edst);
        c.last_pair(&I, &J);
        height[step] = c.ht(I);
        step = c.last_merge(m1+step, m2+step);
        E[step] = e;
    }

    //compute the return values for merge and order

    E[n-1] = 0.0;
    for (i=0; i<n-1; i++) {
        merge[i] = m1[i];
        merge[i+n-1] = m2[i]; //n-1 by 2, stored by columns
    }

    c.order(order,1);

    Free(E);
    Free(m1);
    Free(m2);
    free_matrix(dst, n, n);
    free_matrix(Edst, n, n);
    return;
}

void lower2square(double **dst, double *diss, int n) {
    // convert lower.tri in vector form to square matrix
    // an R distance vector will be in column order by default

    int i, j, k=0;
    for (i=0; i<n; i++)
        dst[i][i] = 0.0;
    k = 0;
    for (j=0; j<n; j++)
        for (i=j+1; i<n; i++) {
            dst[i][j] = dst[j][i] = diss[k];
            k++;
        }
    return;
}

} // end extern "C"
