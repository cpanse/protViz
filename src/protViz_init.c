#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void __findNN(void *, void *, void *, void *, void *);
extern void __findNN_(void *, void *, void *, void *, void *);
extern void _computeFragmentIons(void *, void *, void *, void *);
extern void computeFragmentIons(void *, void *, void *, void *, void *);
extern void computeFragmentIonsModification(void *, void *, void *, void *, void *, void *, void *);
extern void computeParentIonMass(void *, void *, void *);
extern void computeParentIonMass2(void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _protViz_fcat_FASTA(SEXP);
extern SEXP _protViz_tryptic_digest_FASTA(SEXP);
extern SEXP aa2mass_main(SEXP, SEXP, SEXP);
extern SEXP deisotoper_main(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"__findNN",                        (DL_FUNC) &__findNN,                        5},
  {"__findNN_",                       (DL_FUNC) &__findNN_,                       5},
  {"_computeFragmentIons",            (DL_FUNC) &_computeFragmentIons,            4},
  {"computeFragmentIons",             (DL_FUNC) &computeFragmentIons,             5},
  {"computeFragmentIonsModification", (DL_FUNC) &computeFragmentIonsModification, 7},
  {"computeParentIonMass",            (DL_FUNC) &computeParentIonMass,            3},
  {"computeParentIonMass2",           (DL_FUNC) &computeParentIonMass2,           5},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"_protViz_fcat_FASTA",           (DL_FUNC) &_protViz_fcat_FASTA,           1},
  {"_protViz_tryptic_digest_FASTA", (DL_FUNC) &_protViz_tryptic_digest_FASTA, 1},
  {"aa2mass_main",                  (DL_FUNC) &aa2mass_main,                  3},
  {"deisotoper_main",               (DL_FUNC) &deisotoper_main,               5},
  {NULL, NULL, 0}
};

void R_init_protViz(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}