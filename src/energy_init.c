#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* declarations to register native routines in this package */ 

/* .C calls */
extern void dCOV(void *, void *, void *, void *, void *, void *, void *);
extern void dCOVtest(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Emin_hclust(void *, void *, void *, void *, void *);
extern void indepE(void *, void *, void *, void *, void *);
extern void indepEtest(void *, void *, void *, void *, void *, void *, void *);
extern void ksampleEtest(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void poisMstat(void *, void *, void *);

/* .Call calls */
extern SEXP _energy_D_center(SEXP);
extern SEXP _energy_dcovU_stats(SEXP, SEXP);
extern SEXP _energy_mvnEstat(SEXP);
extern SEXP _energy_partial_dcor(SEXP, SEXP, SEXP);
extern SEXP _energy_partial_dcov(SEXP, SEXP, SEXP);
extern SEXP _energy_projection(SEXP, SEXP);
extern SEXP _energy_U_center(SEXP);
extern SEXP _energy_U_product(SEXP, SEXP);
extern SEXP _energy_dcov_UV(SEXP, SEXP, SEXP);
extern SEXP _energy_sum_paired_dist(SEXP, SEXP);
extern SEXP _energy_sum_dist3(SEXP, SEXP);
extern SEXP _energy_rowSums_dist(SEXP);
extern SEXP _energy_Btree_sum(SEXP, SEXP);
extern SEXP _energy_get_dcov_sums(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"dCOV",         (DL_FUNC) &dCOV,         7},
  {"dCOVtest",     (DL_FUNC) &dCOVtest,     8},
  {"Emin_hclust",  (DL_FUNC) &Emin_hclust,  5},
  {"indepE",       (DL_FUNC) &indepE,       5},
  {"indepEtest",   (DL_FUNC) &indepEtest,   7},
  {"ksampleEtest", (DL_FUNC) &ksampleEtest, 9},
  {"poisMstat",    (DL_FUNC) &poisMstat,    3},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"_energy_D_center",       (DL_FUNC) &_energy_D_center,      1},
  {"_energy_dcovU_stats",    (DL_FUNC) &_energy_dcovU_stats,   2},
  {"_energy_mvnEstat",       (DL_FUNC) &_energy_mvnEstat,      1},
  {"_energy_partial_dcor",   (DL_FUNC) &_energy_partial_dcor,  3},
  {"_energy_partial_dcov",   (DL_FUNC) &_energy_partial_dcov,  3},
  {"_energy_projection",     (DL_FUNC) &_energy_projection,    2},
  {"_energy_U_center",       (DL_FUNC) &_energy_U_center,      1},
  {"_energy_U_product",      (DL_FUNC) &_energy_U_product,     2},
  {"_energy_dcov_UV",        (DL_FUNC) &_energy_dcov_UV,       3},
  {"_energy_sum_paired_dist",(DL_FUNC) &_energy_sum_paired_dist,  2},
  {"_energy_sum_dist3",      (DL_FUNC) &_energy_sum_dist3,     2},
  {"_energy_rowSums_dist",   (DL_FUNC) &_energy_rowSums_dist,  1},
  {"_energy_Btree_sum",      (DL_FUNC) &_energy_Btree_sum,     2},
  {"_energy_get_dcov_sums",  (DL_FUNC) &_energy_get_dcov_sums, 3},
  {NULL, NULL, 0}
};

void R_init_energy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
