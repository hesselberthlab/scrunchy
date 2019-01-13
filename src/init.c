#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _scrunchy_compute_snn_impl(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_scrunchy_compute_snn_impl", (DL_FUNC) &_scrunchy_compute_snn_impl, 2},
    {NULL, NULL, 0}
};

void R_init_scrunchy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
