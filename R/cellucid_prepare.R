# The help page for cellucid_prepare() is hand-written in
# man/cellucid_prepare.Rd and is the only source for it. Nothing in this
# package is generated: NAMESPACE is hand-written too, because
# useDynLib(cellucid, .registration = TRUE, .fixes = "C_") has no generator
# input anywhere under R/. Edit the .Rd file directly when this signature
# changes; tests/testthat/test-current-contract.R rebuilds the documented
# \usage block into a function and compares it to these formals, so the two
# cannot drift apart.
cellucid_prepare <- function(
    out_path,
    format,
    data = NULL,
    embeddings = NULL,
    dataset_id = NULL,
    dataset_name = NULL,
    obs = NULL,
    var = NULL,
    gene_expression = NULL,
    connectivities = NULL,
    vector_fields = NULL,
    modalities = NULL,
    gene_identifiers = NULL,
    dataset_description = NULL,
    citation_text = NULL,
    citation_url = NULL,
    source_name = NULL,
    source_url = NULL,
    embedding_defaults = NULL,
    latent_key = NULL,
    obs_keys = NULL,
    var_gene_id_column = NULL,
    vector_field_default = NULL,
    normalize_embeddings = NULL,
    centroid_min_points = NULL,
    modality = NULL,
    zarr_format = NULL,
    force = FALSE
) {
  # An argument that was never supplied and an argument supplied with the wrong
  # value are two different mistakes, and only the message tells the caller
  # which one they made. Every required argument is therefore reported here,
  # before any value is inspected. Each `missing()` sits directly in an `if`
  # condition in this body, which is the one place R guarantees it reads the
  # argument's own supplied-ness. Without `data`, the embeddings and the
  # dataset's identity are required too; with it, they default to what the
  # object holds (.plan_dataset() names what is still missing).
  absent_required <- character()
  if (missing(out_path)) {
    absent_required <- c(absent_required, "out_path")
  }
  if (missing(format)) {
    absent_required <- c(absent_required, "format")
  }
  if (is.null(data)) {
    if (is.null(embeddings)) {
      absent_required <- c(absent_required, "embeddings")
    }
    if (is.null(dataset_id)) {
      absent_required <- c(absent_required, "dataset_id")
    }
    if (is.null(dataset_name)) {
      absent_required <- c(absent_required, "dataset_name")
    }
  }
  .require_supplied_arguments(absent_required)

  format <- .validate_format(format)
  zarr_format <- .validate_zarr_format(zarr_format)
  force <- .validate_flag(force, "force")
  target <- .resolve_output_path(out_path, format = format, force = force)
  target$zarr_format <- zarr_format

  components <- .dataset_components(data, list(
    embeddings = embeddings,
    dataset_id = dataset_id,
    dataset_name = dataset_name,
    obs = obs,
    var = var,
    gene_expression = gene_expression,
    connectivities = connectivities,
    vector_fields = vector_fields,
    modalities = modalities,
    gene_identifiers = gene_identifiers,
    dataset_description = dataset_description,
    citation_text = citation_text,
    citation_url = citation_url,
    source_name = source_name,
    source_url = source_url,
    embedding_defaults = embedding_defaults,
    latent_key = latent_key,
    obs_keys = obs_keys,
    var_gene_id_column = var_gene_id_column,
    vector_field_default = vector_field_default,
    modality = modality,
    normalize_embeddings = normalize_embeddings,
    centroid_min_points = centroid_min_points
  ))
  plan <- .plan_dataset(components)
  .write_output(target, write = plan$write)
  invisible(target$path)
}
