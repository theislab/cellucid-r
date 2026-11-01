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
SEXP cellucid_gene_major_block(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers, SEXP cells,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype, SEXP threads);
SEXP cellucid_renumbered_block(SEXP i, SEXP p, SEXP x, SEXP new_row, SEXP column_order, SEXP pointers,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype, SEXP part, SEXP padded);
SEXP cellucid_raw_slice(SEXP bytes, SEXP offset, SEXP count, SEXP total);
SEXP cellucid_layouts_start(SEXP jobs, SEXP paths, SEXP threads, SEXP park);
SEXP cellucid_layouts_unpark(SEXP pointer);
SEXP cellucid_layouts_stop(SEXP pointer);
SEXP cellucid_layout_state(SEXP pointer, SEXP index);
SEXP cellucid_layout_slice(SEXP pointer, SEXP index, SEXP part, SEXP offset, SEXP count, SEXP total);
SEXP cellucid_layout_block(SEXP pointer, SEXP index, SEXP first, SEXP last);
SEXP cellucid_cells_frame(SEXP header, SEXP cells, SEXP allow_empty);
SEXP cellucid_backed_tally(SEXP n_rows);
SEXP cellucid_backed_count(SEXP tally, SEXP index, SEXP pointers, SEXP values, SEXP by_rows, SEXP rows, SEXP col0);
SEXP cellucid_backed_counted(SEXP tally);
SEXP cellucid_backed_prefix_tally(SEXP n_genes, SEXP levels);
SEXP cellucid_backed_prefix_count(SEXP tally, SEXP index, SEXP pointers, SEXP values, SEXP by_rows, SEXP genes, SEXP col0,
                                  SEXP new_row);
SEXP cellucid_backed_prefix_counted(SEXP tally);
SEXP cellucid_backed_target(SEXP pointers, SEXP low, SEXP high, SEXP n_cells, SEXP index_size, SEXP dtype, SEXP paths);
SEXP cellucid_backed_target_genes(SEXP target);
SEXP cellucid_backed_place(SEXP target, SEXP first_gene, SEXP last_gene, SEXP index, SEXP pointers, SEXP values,
                           SEXP by_rows, SEXP genes, SEXP col0, SEXP new_row, SEXP threads);
SEXP cellucid_backed_finish(SEXP target, SEXP first_gene, SEXP last_gene, SEXP levels, SEXP threads);
SEXP cellucid_backed_slice(SEXP target, SEXP part, SEXP first, SEXP count, SEXP total);
SEXP cellucid_backed_release(SEXP target);
SEXP cellucid_subsample_rank(SEXP seed, SEXP positions);
SEXP cellucid_prefix_counts(SEXP p, SEXP i, SEXP levels, SEXP n_genes, SEXP cells_by_column, SEXP map,
                            SEXP cells, SEXP new_row);
SEXP cellucid_subsample_uniform(SEXP n_obs, SEXP k, SEXP seed);
SEXP cellucid_subsample_blocks(SEXP n_obs, SEXP k, SEXP seed, SEXP block_size);
SEXP cellucid_subsample_split(SEXP strata, SEXP n_strata, SEXP k, SEXP seed, SEXP allocation, SEXP minimum);
SEXP cellucid_subsample_boxes(SEXP coordinates, SEXP dimensions, SEXP grid_size);
SEXP cellucid_subsample_digest(SEXP n_obs, SEXP kept);
SEXP cellucid_subsample_encode(SEXP kept);
SEXP cellucid_subsample_decode(SEXP payload, SEXP n_obs);
SEXP cellucid_shortest_decimal(SEXP x);
SEXP cellucid_file_version(SEXP path);
SEXP cellucid_file_range(SEXP path, SEXP offset, SEXP length);
SEXP cellucid_deflate_chunks(SEXP chunks, SEXP level, SEXP threads);
SEXP cellucid_process_alive(SEXP pid);
SEXP cellucid_physical_memory(void);

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
    {"cellucid_gene_major_block", (DL_FUNC)&cellucid_gene_major_block, 12},
    {"cellucid_renumbered_block", (DL_FUNC)&cellucid_renumbered_block, 12},
    {"cellucid_raw_slice", (DL_FUNC)&cellucid_raw_slice, 4},
    {"cellucid_layouts_start", (DL_FUNC)&cellucid_layouts_start, 4},
    {"cellucid_layouts_unpark", (DL_FUNC)&cellucid_layouts_unpark, 1},
    {"cellucid_layouts_stop", (DL_FUNC)&cellucid_layouts_stop, 1},
    {"cellucid_layout_state", (DL_FUNC)&cellucid_layout_state, 2},
    {"cellucid_layout_slice", (DL_FUNC)&cellucid_layout_slice, 6},
    {"cellucid_layout_block", (DL_FUNC)&cellucid_layout_block, 4},
    {"cellucid_cells_frame", (DL_FUNC)&cellucid_cells_frame, 3},
    {"cellucid_backed_tally", (DL_FUNC)&cellucid_backed_tally, 1},
    {"cellucid_backed_count", (DL_FUNC)&cellucid_backed_count, 7},
    {"cellucid_backed_counted", (DL_FUNC)&cellucid_backed_counted, 1},
    {"cellucid_backed_prefix_tally", (DL_FUNC)&cellucid_backed_prefix_tally, 2},
    {"cellucid_backed_prefix_count", (DL_FUNC)&cellucid_backed_prefix_count, 8},
    {"cellucid_backed_prefix_counted", (DL_FUNC)&cellucid_backed_prefix_counted, 1},
    {"cellucid_backed_target", (DL_FUNC)&cellucid_backed_target, 7},
    {"cellucid_backed_target_genes", (DL_FUNC)&cellucid_backed_target_genes, 1},
    {"cellucid_backed_place", (DL_FUNC)&cellucid_backed_place, 11},
    {"cellucid_backed_finish", (DL_FUNC)&cellucid_backed_finish, 5},
    {"cellucid_backed_slice", (DL_FUNC)&cellucid_backed_slice, 5},
    {"cellucid_backed_release", (DL_FUNC)&cellucid_backed_release, 1},
    {"cellucid_subsample_rank", (DL_FUNC)&cellucid_subsample_rank, 2},
    {"cellucid_subsample_uniform", (DL_FUNC)&cellucid_subsample_uniform, 3},
    {"cellucid_subsample_blocks", (DL_FUNC)&cellucid_subsample_blocks, 4},
    {"cellucid_subsample_split", (DL_FUNC)&cellucid_subsample_split, 6},
    {"cellucid_subsample_boxes", (DL_FUNC)&cellucid_subsample_boxes, 3},
    {"cellucid_subsample_digest", (DL_FUNC)&cellucid_subsample_digest, 2},
    {"cellucid_subsample_encode", (DL_FUNC)&cellucid_subsample_encode, 1},
    {"cellucid_subsample_decode", (DL_FUNC)&cellucid_subsample_decode, 2},
    {"cellucid_shortest_decimal", (DL_FUNC)&cellucid_shortest_decimal, 1},
    {"cellucid_prefix_counts", (DL_FUNC)&cellucid_prefix_counts, 8},
    {"cellucid_file_version", (DL_FUNC)&cellucid_file_version, 1},
    {"cellucid_file_range", (DL_FUNC)&cellucid_file_range, 3},
    {"cellucid_deflate_chunks", (DL_FUNC)&cellucid_deflate_chunks, 3},
    {"cellucid_process_alive", (DL_FUNC)&cellucid_process_alive, 1},
    {"cellucid_physical_memory", (DL_FUNC)&cellucid_physical_memory, 0},
    {NULL, NULL, 0}};

void R_init_cellucid(DllInfo *dll) {
  R_registerRoutines(dll, NULL, cellucid_call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
