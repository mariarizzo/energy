/*
   ECl.cc:   energy package

   Author:  Maria Rizzo <mrizzo @ bgnet.bgsu.edu>
   Created: 12 Dec 2002
   Revised: 4 Jan 2004 for R-1.8.1 energy package
   Revised: 28 Jan 2004

*/

#include "ECl.h"
#include <R.h>

extern "C" {

//implementation of Cl and ECl classes

Cl::~Cl() {  //destructor
    int i;
    if (isinit==1) {
        Free(size);
        Free(step);
        Free(height);
        Free(w);
        for (i=0; i<n; i++)
            Free(clus[i]);
        Free(clus);
    }
}

int Cl::init(int m) {
    int i, j;

    if (isinit==1 && m!=n) error("is initialized");
    n = m;
    nclus = m;
    if (isinit==0) {
        size =  Calloc(n, int);
        step = Calloc(n, int);
        height = Calloc(n, double);
        w = Calloc(n, int);
        clus = Calloc(n, int *);
        for (i = 0; i < n; i++)
            clus[i] = Calloc(n, int);
        }
    for (i=0; i<n; i++) {
        size[i] = 1;
        step[i] = -i-1;
        height[i] = -1;
        for (j=0; j<n; j++)
            clus[i][0] = i;
    }
    r1 = r2 = c1 = c2 = n;
    it = 0;
    isinit = 1;
    return nclus;
}

int Cl::init(int m, int *m1, int *m2, int k) {
    //initialize to k groups from merge order, indexed from 1
    int i, I, J;

    init(m);
    if (k < n && k > 0) {
        I = -m1[0]-1;
        J = -m2[0]-1;
        combine(I, J);
        w[0] = J;
        w[1] = I;
        i = 1;
        while (nclus > k) {
            I = (m1[i]<0)? -m1[i]-1 : w[m1[i]];
            J = (m2[i]<0)? -m2[i]-1 : w[m2[i]];
            combine(I, J);
            i++;
            w[i] = I;
        }
    }
    nclus = clusters();
    return nclus;
}

int Cl::init(int m, int *G, int base) {
    //initialize cluster with group membership vector G
    int g, i;

    init(m);
    if (base > 0)
        for (i=0; i<n; i++) G[i]-=base;
    for (i=0; i<n; i++) w[i] = 0;
    for (i=0; i<n; i++) {
        g = G[i];
        clus[g][w[g]] = i;
        w[g]++;
    }
    for (i=0; i<n; i++) {
        size[i] = w[i];
        step[i] = 0;
        height[i] = -1;
    }
    r1=c1=r2=c2=12;
    pstep1=pstep2=psize1=psize2=0;
    nclus = clusters();
    return nclus;
}


int Cl::clusters(int *cl) {
    //cl is indices of non-empty clusters
    //count non-empty clusters and check value of nclus and n
    int i, k=0, m=0;
    for (i=0; i<n; i++)
        if (size[i]>0) {
            cl[k++] = i;
            m+= size[i];
        }
    if (k!=nclus) error("nclus error");
    if (m!=n) error("total size error");
    return nclus;
}

int Cl::clusters() {
    //count the number of non-empty clusters, and set nclus
    int i, k=0;
    for (i=0; i<n; i++)
        if (size[i]>0) k++;
    if (k>n || k<1) error("nclus error");
    nclus = k;
    return nclus;
}

int Cl::combine(int I, int J) {
    //merge Jth row into Ith row
    //w is preserved
    int j, m;

    if (I==J) error("c:I==J");
    if (I<0 || J<0 || I>=n || J>=n) error("c:I,J error");
    if (size[I]<=0 || size[J]<=0) error("c:empty cluster");
    if (nclus < 2) error("c:1 cluster");

    m = size[I];
    for (j=0; j<size[J]; j++)
        clus[I][m+j] = clus[J][j];
    psize1 = size[I];
    psize2 = size[J];
    size[I] += size[J];
    size[J] = 0;
    nclus--;
    pstep1 = step[I];
    step[I] = n - nclus;
    pstep2 = step[J];
    r1 = I; r2 = J; //cannot be undetached
    c1 = c2 = n;
    it++;
    return 1;
}

int Cl::groups(int *g, int base) {
    //current group membership:  x[i] in g[i],
    //return number of groups
    int i, j, k=0;
    for (i=0; i<n; i++)
        if (size[i] > 0) {
            for (j=0; j<size[i]; j++)
                g[clus[i][j]] = k;
            k++;
        }
    if (base > 0)
        for (i=0; i<n; i++)
            g[i]+= base;
    return nclus;
}

int Cl::last_merge(int *s1, int *s2) {
    //[s1,s2] for merge matrix, hierarchical clustering
    //return value is the step
    *s1 = pstep1;
    *s2 = pstep2;
    return n-nclus;
}

int Cl::last_pair(int *i, int *j) {
    //last pair combined
    //return value is the step
    *i = r1;
    *j = r2;
    return n-nclus;
}

int Cl::order(int *o, int base) {
    //permutation of indices such that branches in a
    //hierarchical tree will not cross
    int i, j, k=0;
    for (i=0; i<n; i++)
        if (size[i] > 0)
            for (j=0; j<size[i]; j++)
                o[k++] = clus[i][j];
    if (base > 0)
        for (i=0; i<n; i++) o[i]+= base;
    if (k > n) return -1;
    return 0;
}

int Cl::proximity(int **p) {
    //p[i][j] is 1 if (i,j) in same cluster
    //p[i][j] is 0 if (i,j) in different clusters
    int a, b, i, j, k;
    for (i=0; i<n; i++) {
        p[i][i] = 1;
        for (j=i+1; j<n; j++)
            p[i][j] = p[j][i] = 0;
    }
    for (i=0; i<n; i++)
        for (j=0; j<size[i]; j++)
            for (k=0; k<j; k++)
            {
                a = clus[i][j];
                b = clus[i][k];
                p[a][b] = p[b][a] = 1;
            }
    return nclus;
}

////////////////////////////////////////////////////////////////////////////

double ECl::cldst(int k1, int k2, double **dst)
{
    // return the E-distance between clusters k1, k2
    int     a, b, i, j, m, n;
    double  e, sumxy, sumxx, sumyy, h;

    m = size[k1];
    n = size[k2];
    if (m==0 || n==0) return 0.0;
    if (k1==k2) return 0.0;

    sumxy = 0.0;
    for (i=0; i<m; i++) {
        a = clus[k1][i];
        for (j=0; j<n; j++) {
            b = clus[k2][j];
            sumxy += dst[a][b];
        }
    }
    sumxx = 0.0;
    for (i=0; i<m; i++) {
        a = clus[k1][i];
        for (j=0; j<i; j++) {
            b = clus[k1][j];
            sumxx += dst[a][b];
        }
    }
    sumxx *= 2.0;
    sumyy = 0.0;
    for (i=0; i<n; i++) {
        a = clus[k2][i];
        for (j=0; j<i; j++) {
            b = clus[k2][j];
            sumyy += dst[a][b];
        }
    }
    sumyy *= 2.0;
    h = (double)(2*m*n) / (double)(m+n);
    e = h*( 2.0*sumxy/(double)(m*n) -
             sumxx/(double)(m*m) - sumyy/(double)(n*n) );
    return e;
}

double ECl::calc_E(double **dst) {
    //sum cluster distances over all non-empty clusters = E
    //directly from the dissimilarity matrix
    //faster to use Edst stored matrix

    int    i, j;
    double sum=0.0;

    for (i=next_cl(-1); i<n;  ) {
        for (j=next_cl(i); j<n;  ) {
            sum += cldst(i, j, dst);
            j=next_cl(j);
        }
        i=next_cl(i);
    }
    pE = E;
    E = sum;
    return sum;
}

double ECl::sum_Edst(double **Ed) {
    //sum E-distances over all non-empty clusters
    int    i, j;
    double sum=0.0;

    for (i=next_cl(-1); i<n;  ) {
        for (j=next_cl(i); j<n;  ) {
            sum += Ed[i][j];
            j=next_cl(j);
        }
        i=next_cl(i);
    }
    pE = E;
    E = sum;
    return sum;
}

double ECl::update_Edst(int cl1, int cl2, double **dst, double **Ed) {
    //update E-distance array after merging clusters cl1, cl2
    int k;

    for (k=0; k<n; k++) {
        Ed[cl1][k] = Ed[k][cl1] = cldst(cl1, k, dst);
        Ed[cl2][k] = Ed[k][cl2] = cldst(cl1, k, dst);
    }
    return sum_Edst(Ed);
}

double ECl::update_Edst(double **dst, double **Ed) {
    //recalculate all in E-distance array
    int i, j;

    for (i=0; i<n; i++) {
        Ed[i][i] = 0.0;
        for (j=0; j<n; j++)
            Ed[i][j] = Ed[j][i] = cldst(i, j, dst);
    }
    return sum_Edst(Ed);
}

double ECl::init_Edst(double **dst, double **Ed) {
    int    i, j;
    E=0.0;
    for (i=0; i<n; i++) {
        Ed[i][i] = 0.0;
        for (j=i+1; j<n; j++) {
            Ed[i][j] = Ed[j][i] = 2.0*dst[i][j];
            E += Ed[i][j];
        }
    }
    pE = E;
    return E;
}

double ECl::find_minEdst(double **Ed, int *imin, int *jmin) {
    //find min E distance over all non-empty clusters

    int    i, j, I, J;
    double m;

    I=i=next_cl(-1);
    J=j=next_cl(i);
    m=Ed[i][j];
    for (i=next_cl(-1); i<n;  ) {
        for (j=next_cl(i); j<n;  ) {
            if (Ed[i][j] < m) {
                m=Ed[i][j];
                I=i;
                J=j;
            }
            j=next_cl(j);
        }
        i=next_cl(i);
    }
    (*imin)=I;
    (*jmin)=J;
    return m;
}


double ECl::merge_minEdst(double **dst, double **Ed) {
    //merge the pair of clusters that have the minimum
    //pairwise E-cluster distance
    int    d, p, I, J;
    double hI, hJ;

    clusters(w);
    if (nclus == 2) {
        I=w[0]; J=w[1];
        if (height[I] > height[J]) {
            I=w[1]; J=w[0];
        }
        height[I]=Ed[I][J];
        combine(I, J);
        update_Edst(I, J, dst, Ed);
        return 0.0;
    }

    if (nclus == 1) error("last cluster");
    if (nclus < 1) error("nclus<1");

    I=J=-1;
    find_minEdst(Ed, &I, &J);

    if (I>=0) {
        if (J < I) {
            p=I; I=J; J=p;
        }
        hI = hJ = 0.0;
        if (step[I] > 0) hI = height[I];
        if (step[J] > 0) hJ = height[J];
        if (hJ < hI) {
            p=I; I=J; J=p;
        }

        height[I] = Ed[I][J];
        d=combine(I,J);
        if(!d) error("merge_best_pair error");
        pE = E;
        E = update_Edst(I, J, dst, Ed);
    }
    return E;
}


} //end extern "C"
