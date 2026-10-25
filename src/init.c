/*
 * The registration of every native routine of the package. NAMESPACE loads
 * them with useDynLib(cellucid, .registration = TRUE, .fixes = "C_"), and
 * R/native.R is their one R-side caller.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP cellucid_crc32(SEXP bytes, SEXP previous);
SEXP cellucid_random_bytes(SEXP n);
SEXP cellucid_sha256(SEXP bytes);
SEXP cellucid_extent(SEXP values, SEXP check_whole);
SEXP cellucid_row_pointers(SEXP i, SEXP map, SEXP n_rows);
SEXP cellucid_kept_extent(SEXP x, SEXP i, SEXP map, SEXP check_whole);
SEXP cellucid_latent_outlier_quantiles(SEXP latent, SEXP codes,
                                       SEXP n_categories, SEXP min_points);
SEXP cellucid_inflate(SEXP bytes, SEXP gzip, SEXP expected);
SEXP cellucid_digest_cell_set(SEXP indices);
SEXP cellucid_decode_cell_set(SEXP payload, SEXP cell_count, SEXP reference);
SEXP cellucid_decode_user_codes(SEXP payload, SEXP expected_length, SEXP codes_type);
SEXP cellucid_encode_values(SEXP values, SEXP first, SEXP count, SEXP chunk_length, SEXP dtype);
SEXP cellucid_gene_major_block(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype);
SEXP cellucid_raw_slice(SEXP bytes, SEXP offset, SEXP count, SEXP total);
SEXP cellucid_cells_frame(SEXP header, SEXP cells, SEXP allow_empty);

static const R_CallMethodDef cellucid_call_methods[] = {
    {"cellucid_crc32", (DL_FUNC)&cellucid_crc32, 2},
    {"cellucid_random_bytes", (DL_FUNC)&cellucid_random_bytes, 1},
    {"cellucid_sha256", (DL_FUNC)&cellucid_sha256, 1},
    {"cellucid_extent", (DL_FUNC)&cellucid_extent, 2},
    {"cellucid_row_pointers", (DL_FUNC)&cellucid_row_pointers, 3},
    {"cellucid_kept_extent", (DL_FUNC)&cellucid_kept_extent, 4},
    {"cellucid_latent_outlier_quantiles", (DL_FUNC)&cellucid_latent_outlier_quantiles, 4},
    {"cellucid_inflate", (DL_FUNC)&cellucid_inflate, 3},
    {"cellucid_digest_cell_set", (DL_FUNC)&cellucid_digest_cell_set, 1},
    {"cellucid_decode_cell_set", (DL_FUNC)&cellucid_decode_cell_set, 3},
    {"cellucid_decode_user_codes", (DL_FUNC)&cellucid_decode_user_codes, 3},
    {"cellucid_encode_values", (DL_FUNC)&cellucid_encode_values, 5},
    {"cellucid_gene_major_block", (DL_FUNC)&cellucid_gene_major_block, 10},
    {"cellucid_raw_slice", (DL_FUNC)&cellucid_raw_slice, 4},
    {"cellucid_cells_frame", (DL_FUNC)&cellucid_cells_frame, 3},
    {NULL, NULL, 0}};

void R_init_cellucid(DllInfo *dll) {
  R_registerRoutines(dll, NULL, cellucid_call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
