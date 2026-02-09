#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP conc(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"conc", (DL_FUNC) &conc, 4},
  {NULL, NULL, 0}
};

void R_init_RCor(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
