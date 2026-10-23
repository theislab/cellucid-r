/*
 * The registration of every native routine of the package. NAMESPACE loads
 * them with useDynLib(cellucid, .registration = TRUE, .fixes = "C_"), and
 * R/native.R is their one R-side caller.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP cellucid_crc32(SEXP bytes, SEXP previous);
SEXP cellucid_latent_outlier_quantiles(SEXP latent, SEXP codes,
                                       SEXP n_categories, SEXP min_points);

static const R_CallMethodDef cellucid_call_methods[] = {
    {"cellucid_crc32", (DL_FUNC)&cellucid_crc32, 2},
    {"cellucid_latent_outlier_quantiles", (DL_FUNC)&cellucid_latent_outlier_quantiles, 4},
    {NULL, NULL, 0}};

void R_init_cellucid(DllInfo *dll) {
  R_registerRoutines(dll, NULL, cellucid_call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
