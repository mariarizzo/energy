/*
   ECl.h:   energy package

   Author:  Maria Rizzo <mrizzo @ bgnet.bgsu.edu>
   Created: 12 Dec 2002
   Revised: 4 Jan 2004 for R-1.8.1 energy package
   Revised: 28 Jan 2004
*/

#define DLLEXPORT

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<float.h>

//declarations for class Cl for hierarchical cluster analysis
//and ECl for e-clustering

#define EPS DBL_EPSILON*20.0
#define ONE 1.0+DBL_EPSILON*20.0

class DLLEXPORT Cl {
    //for cluster analysis
protected:
    int    n;        //number of observations
    int    nclus;    //number of clusters
    int    it;       //number of changes to clusters
    int    pstep1;
    int    pstep2;
    int    psize1;
    int    psize2;
    int    r1;
    int    r2;
    int    c1;
    int    c2;
    int    temp;
    int    isinit;   //is memory allocated for arrays
    int    *size;    //sizes of clusters
    int    *step;    //step when cluster formed
    double *height;  //distance between merging clusters
    int    *w;
    int    **clus;   //indices of observations

public:
    Cl(){isinit=0;};  //no memory is allocated, call init(n)
    ~Cl();
    int    init(int m);
    int    init(int n, int *m1, int *m2, int k);
    int    init(int m, int *G, int base);

    int    dim() {return n;}
    int    len() {return nclus;}
    int    len(int i) {return size[i];}
    int    obs(int i, int j) {return clus[i][j];}
    double ht(int i) {return height[i];}
    int    next_cl(int p) {p++;while(p<n && size[p]<1) p++;return p;}
    int    clusters();
    int    clusters(int *cl);
    int    combine(int I, int J);
    int    last_merge(int *s1, int *s2);
    int    last_pair(int *i, int *j);
    int    groups(int *g, int base);
    int    order(int *o, int base);
    int    proximity(int **p);

};

class DLLEXPORT ECl:public Cl {
    //derived class supports E-clustering algorithm
protected:
    double E;
    double pE;
public:
    double cldst(int cl1, int cl2, double **dst);
    double calc_E(double **dst);
    double sum_Edst(double **Ed);
    double init_Edst(double **dst, double **Ed);
    double update_Edst(double **dst, double **Ed);
    double update_Edst(int cl1, int cl2, double **dst, double **Ed);
    double find_minEdst(double **Ed, int *imin, int *jmin);
    double merge_minEdst(double **dst, double **Ed);
    double merge_best_pair(double **dst, double **Ed);
};

