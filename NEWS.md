# cellucid 0.9.5 <!-- CELLUCID_VERSION -->

Version 0.9.5 is the CRAN submission release.

It is a patch release on 0.9.4 and stays in step with the Python package of
the same version. It breaks with every earlier release: the prepared export
format is removed, and `cellucid_prepare()` writes a Cellucid-ready AnnData
instead, a standard AnnData Zarr store that Python's anndata reads and the
Cellucid web app opens directly. The contract is the documentation's
"Cellucid-ready AnnData specification" page. An export written by 0.9.4 is not
read by the matching web app; prepare the data again with this release.

- `cellucid_prepare(out_path, format, embeddings, dataset_id, dataset_name,
  ...)` takes a required `format` right after `out_path`: `"zarr"` (a `.zarr`
  folder of small files, the way to publish on GitHub: no file over 100 MB, no
  LFS) or `"zarr.zip"` (the same store as one `.zarr.zip` file, for hosts that
  serve single files or a local file). It has no default: a missing `format`
  is reported with the other missing required arguments, naming both choices.
  A value is matched as `match.arg()` matches, anything else is refused with
  the choices, and `out_path` must end in the suffix of `format`.
  `zarr_format = 3` writes Zarr v3 (`zarr.json` documents, `c/...` chunk keys,
  consolidated metadata in the root `zarr.json`); Zarr v2 is the default.
  There is no `"h5ad"`: the package writes no HDF5 and depends on no HDF5
  library; Python's `cellucid.prepare(format = "h5ad")` writes `.h5ad`.
- The store holds `X` gene-major (CSC from a sparse matrix, stored zeros
  dropped) in float32, or float64 when the genes kept hold integer counts
  above 2^24; `obs` and `var` with anndata's own encodings and pandas' category
  code widths (int8 below 127 categories, int16 below 32,767, int32 beyond);
  every embedding whole in `obsm` (float32, unnormalized);
  `obsp/connectivities` as symmetric float32 CSR; and `uns[["cellucid"]]`
  (`schema`, `id`, `name`, `description`, `citation`, `source`, and only the
  launch `defaults` that were set) with `uns[["title"]]`. Chunks hold 2^16
  values for `X` and 2^18 for arrays read whole, deflated at level 6 after a
  byte shuffle on every multi-byte array except the expression values: in
  Zarr v2 numcodecs' `shuffle` filter and zlib; in Zarr v3 only codecs of the
  v3 specification, the `blosc` codec (byte shuffle, zlib level 6, the frames
  c-blosc writes) on shuffled arrays and gzip on the others. A `.zarr.zip` stores chunks first in grid order and
  the consolidated metadata last, with fixed timestamps and permissions. The
  store is written beside its target and moved into place when complete.
- R and Python write identical Zarr v2 and v3 stores from the same input:
  every metadata document (encodings, dtypes, shapes, chunks, codecs, fill
  values) and every value, checked by `tests/testthat/test-cross-writer.R`.
  Every `NaN` and `NA` in a float array is written as numpy's canonical quiet
  NaN on every platform (R's `NaN` has its sign bit set on x86-64).
- New arguments: `out_path`, `format`, `zarr_format`, `citation_text`,
  `citation_url`, `latent_key` (the `embeddings` entry the outlier filter
  measures in) and `normalize_embeddings`. Choose the columns an embedding
  opens on with `embedding_defaults = list("3d" = list(key = "X_pca", columns
  = 1:3))`.
- With `latent_key`, `cellucid_prepare()` stores the latent-space outlier
  quantiles of every factor, character or logical `obs` column in
  `uns[["cellucid"]][["outlierQuantiles"]]` (`latentKey`, `minPoints` =
  `centroid_min_points` or 10, and one float32 array per field, `NaN` where a
  cell has none), so the viewer filters outliers without reading the latent
  space. They are computed in C (`src/latent_outliers.c`), byte for byte what
  Python's `prepare()` stores and what the web app computes, pinned by the
  shared conformance fixture
  (`tests/testthat/fixtures/latent_outlier_conformance.json`).
- Vector fields follow scVelo's naming on any embedding: `vector_fields =
  list(velocity = list(X_umap = v))` writes `obsm[["velocity_umap"]]`, with as
  many columns as the embedding. An embedding whose key does not start with
  `X_`, an array of another width, two fields that would share a key, and a
  key the viewer would read as another field or as an embedding are refused
  before anything is written. `vector_field_default` names the field
  (`"velocity"`) and must be drawn on an embedding the viewer shows when the
  store opens.
- New `cellucid_transition_drift(transition_matrix, embedding,
  normalize_rows)` computes each cell's drift in an embedding's coordinates
  from a cell-cell transition matrix (CellRank's, for example): `T %*% E - E`,
  with the rows of `T` divided by their sums first when `normalize_rows =
  TRUE` (required; a zero row sum is refused). It is Python's
  `compute_transition_drift()` operation for operation, held to the same
  results by a shared conformance fixture
  (`tests/testthat/fixtures/transition_drift_conformance.json`). Store the
  result with `cellucid_prepare(vector_fields = list(T_fwd = list(X_umap =
  drift)))`.
- A sparse neighbour graph's stored zeros are dropped (`Matrix::drop0()`)
  instead of refusing the graph; duplicate coordinates are still refused.
- Removed with the export format: the `latent_space` (name an embedding with
  `latent_key`), `out_dir`, `obs_categorical_dtype`, `var_quantization`,
  `obs_continuous_quantization`, `compression`, `var_gene_name_columns`,
  `source_citation` (now `citation_text` and `citation_url`), `created_at` and
  `centroid_outlier_quantile` arguments, `cellucid_embedding()` (store the
  whole array), and the cross-process export lock and transaction journal.

# cellucid 0.9.4

Version 0.9.4 is the CRAN submission release.

It is a patch release on 0.9.3, which was tagged but never published on CRAN,
stays in step with the Python package of the same version, and writes the same
export format 3: an export written by 0.9.3 does not need to be prepared again.

- `cellucid_prepare(obs_categorical_dtype = "uint32")` writes 32-bit
  categorical codes (`obs/{index}.codes.u32`, missing marker 4,294,967,295)
  for fields with more than 65,535 categories, such as barcodes, byte for byte
  as the Python package does. The error for a field that does not fit the
  declared width names the width that does, for example
  "Pass obs_categorical_dtype='uint32'.".
- `cellucid_prepare()` exports a categorical field where no category reaches
  `centroid_min_points` with quantized outlier quantiles, publishing every
  quantile as missing, instead of stopping partway through the export with
  "it contains no finite values".

# cellucid 0.9.3

Version 0.9.3 is the CRAN submission release.

It carries every change of 0.9.2, which was tagged but never published on
CRAN, and stays in step with the Python package of the same version.

`cellucid_prepare()` writes export format v3, byte for byte the format
`cellucid-python` 0.9.3 writes. Opening a v3 export needs the Cellucid web app
release that reads format 3, and an export written by an earlier version is
regenerated with this one.

- Exports any number of embeddings per dimension. `cellucid_prepare()` takes a
  required named list `embeddings` in place of `X_umap_1d`, `X_umap_2d`, and
  `X_umap_3d`, which are removed. A matrix of 1, 2, or 3 columns is plotted
  whole; the new `cellucid_embedding(values, columns, label, source_key)`
  names the 1-based columns of a wider array (published 0-based), a display
  label, and the key the array came from. `embedding_defaults`, for example
  `c("2d" = "X_umap")`, names the embedding each dimension shows first and is
  required for every dimension with more than one.
- Writes each embedding to `embeddings/{index}_{dimension}d.bin`, ordered by
  dimension and then by id in code-point order, and declares it in
  `dataset_identity.json` version 3 under `embeddings.items` with its label,
  dimension, file, and `source` (`key` and 0-based `columns`). The export root
  holds no payload files any more; `embeddings/` is reconciled against its
  declaration as a sixth surface, beside `obs/`, `var/`, `connectivity/`,
  `vectors/`, and the root itself.
- Binds vector fields to embedding ids: `vector_fields = list(velocity =
  list(X_umap = v2, X_umap_3d = v3))`. Each array is scaled by the factor of
  the embedding it is bound to and written to
  `vectors/{fieldIndex}_{embeddingIndex}.bin`; a field id no longer has to
  follow the `<field>_umap_<n>d` grammar.
- Writes `obs_manifest.json` and `var_manifest.json` as `compact_v2`. A
  categorical entry's centroids are keyed by embedding id
  (`centroidsByEmbedding`).
- Names genes from an integer `var` column by its digits (`7157L` is
  `"7157"`), as the Cellucid viewer and `cellucid serve` read the same column;
  a double or logical column is refused with the reason ("floating-point
  values cannot name genes").
- Computes latent-space outlier quantiles in C with the rule `cellucid-python`
  and the web app share, bit for bit: values rounded to float32 and widened to
  double, centroid sums in cell order, squared distances in column order, the
  share as float32. A shared conformance fixture
  (`tests/testthat/fixtures/latent_outlier_conformance.json`) pins the bytes
  in all three. Quantiles can differ from earlier releases only where two
  distances tie within rounding.
- Accepts repeated gene names. They are made unique the way anndata's
  `var_names_make_unique()` does it -- the second `MATR3` becomes `MATR3-1` --
  with the one rule `cellucid-python` and the web app also apply, and an empty
  name (`NA` or `""`) takes the row's name from `rownames(var)`. What changed
  is reported as a message and recorded in `var_manifest.json`.
  `gene_identifiers` matches the names after they were made unique.
- Adds `var_gene_name_columns`: further `var` columns the viewer offers as
  alternative gene names. Each one is stored as metadata in the new
  `gene_name_sets` list of `var_manifest.json`, which replaces the manifest's
  `var_gene_id_column` key (the argument of that name is unchanged and selects
  the default set); no expression payload is duplicated.

- `citation("cellucid")` now lists the Cellucid preprint first, followed by
  the package manual entry. The preprint is a placeholder until it is published
  on bioRxiv and says so; the entry is generated from the `preferred-citation`
  block that every Cellucid repository's `CITATION.cff` shares.

# cellucid 0.9.2

Version 0.9.2 is the CRAN submission release.

It carries the changes prepared as 0.9.1, which was never published; the number
moved to stay in step with the Python package.

- Exposes `cellucid_prepare()` as the single public preparation entry point.
- Publishes complete prepared generations transactionally and rejects mixed or
  stale output directories.
- Serializes independent R and Python exporters with one persistent exact-target
  lock and recovers ownership after process death without leaking native
  handles.
- Requires exact embedding, observation, gene, vector, and categorical-storage
  inputs before writing any scientific artifact.
- Preserves symmetric positive connectivity weights, exact empty graphs, and
  portable edge-index widths without reinterpreting the source graph.
- Rejects non-finite or non-representable float32 values and duplicate
  scientific identifiers before publication, and says how many and where: a
  refusal counts the offending values and names the first few positions.
  `cellucid-python` reports the same counts; each writer prints positions the
  way its own caller indexes them, and names the categories its own language
  distinguishes.
- Applies the float32 range check wherever a value can enter, not only at the
  write. A quantized field below the smallest subnormal was published as zero,
  and an embedding or `latent_space` coordinate was normalized into range before
  the writer saw it; both ends are now one predicate, checked at the source, as
  in `cellucid-python`.
- Names every payload by its integer position on its axis, so no exported
  filename carries dataset content and every dataset has the same directory
  listing: `var/0.values.f32`, `obs/1.codes.u8`, `vectors/0_2d.bin`, and the
  fixed neutral names under `connectivity/`. Each manifest entry declares its
  own index as its first element -- `[index, name]` and
  `[index, name, minValue, maxValue]` for `fields`, `[index, key]` and
  `[index, key, minValue, maxValue]` for `_continuousFields`, and
  `[index, key, categories, dtype, sentinel, centroids]` with the two outlier
  bounds appended for `_categoricalFields` -- and `pathPattern`,
  `codesPathPattern`, and `outlierPathPattern` substitute `{index}`. The writer
  proves that each axis directory uses exactly the indices `0` through `N-1`
  once each, and that the directory holds exactly the payloads its manifest
  declares, before the generation is published. `obs/` is written by both
  manifest arrays, so those two arrays share one index space. The `cellucid`
  Python package writes the identical layout.
- Reconciles the export root against the generation as well, so the
  `points_<n>d.bin` coordinate payloads are covered by the same guarantee as
  every other payload. They are declared in `dataset_identity.json` rather than
  in an axis manifest, and their declared name and their written name are two
  independent spellings of the `compression` setting, so an export that
  published coordinates under a name the viewer is not told to fetch used to
  succeed and then fail in the browser with no points. The root must now hold
  exactly the manifests, the axis directories, and the declared point payloads
  the generation created, and nothing else.
- Writes `_obsSchemas` as a JSON object in every generation, including one that
  carries no observation field at all. `jsonlite` renders an empty *unnamed*
  list as `[]`, so an export prepared from an `obs` with no columns, or with
  `obs_keys = character(0)`, published `"_obsSchemas":[]`; the viewer calls
  `requireRecord()` on that field and refused the dataset with "expected an
  object", while the `cellucid` Python package wrote `{}` from the same input
  and loaded. The manifest validator now asserts the kind rather than only the
  contents, so the two writers cannot part company here again.
- Reads every manifest back out of the staging directory before publishing and
  requires it to parse to exactly the payload that was validated. Everything
  else in the package validates the payload in memory, which proves nothing
  about the file a reader opens: `jsonlite::toJSON()` writes a non-finite double
  as the JSON *string* `"NaN"`, and an empty unnamed list as an array, so a
  manifest could be well-formed, publish cleanly, and still say something other
  than what was checked. A file that is not valid UTF-8, does not parse, or
  disagrees with the payload on any node's kind, keys, length, or value now
  fails the export and names the node, instead of reaching a browser.
- Drops every filename rule from gene names, `obs_keys`, and vector field ids,
  because none of them names a file any more: no portable-ASCII restriction, no
  case-insensitive collision rule, and no Windows device-name rule. `HLA-DRB1/2`
  and `CON` are exported and recorded exactly as supplied. What each identifier
  must still be is what it is for -- a non-empty string, distinct within its
  axis, and text the viewer can draw exactly as stored. `dataset_id` names the
  export directory, so it alone remains a portable filename component.
- Accepts any vector field key matching `<field>_umap_<1|2|3>d`, with no
  character restriction on `<field>`, matching the `cellucid` Python package
  exactly. Vector fields are emitted in code-point order of their ids on both
  writers, so the same input receives the same payload index in either
  language.
- Requires every string the viewer prints verbatim to read on screen as the
  value it stores. A string category label, an exported gene name, an exported
  `obs` key, a vector field id, `dataset_name`, `dataset_description`,
  `source_name`, `source_url`, and `source_citation` are rejected when they
  carry a control character, one of the zero-width
  characters `U+200B`, `U+2060`, or `U+FEFF`, or leading or trailing whitespace
  of any kind including `U+00A0` NO-BREAK SPACE. An empty category label is
  rejected, as are two labels in one field that a whitespace-collapsing
  renderer draws identically. Nothing is trimmed: trimming would rewrite an
  annotation the caller never asked to change, and would merge `"Liver "` into
  a separate `"Liver"` category and move cells between them. The message names
  every offending label in the field at once and gives the one-line repair. The
  `cellucid` Python package enforces the identical rule.
- Extends the `dataset_name` check, which previously missed `U+0080`-`U+009F`
  control characters and non-ASCII whitespace such as `U+00A0`, to the same
  shared rule.
- Derives quantized continuous payloads and bounds from the viewer's exact
  float32 value domain.
- Publishes a constant continuous field instead of refusing the export. A gene
  expressed at one level in every exported cell -- very often zero, once an
  atlas is subset to one lineage -- and an `obs` column a subset flattened are
  both ordinary data, and quantization now has a named case for them:
  `minValue == maxValue` with every code `0`, which the viewer decodes back to
  the exact constant rather than to an approximation of it. `.quantize_continuous()`
  takes an explicit branch for it and never derives a scale, so nothing divides
  by `maxValue - minValue`. Native-double variation finer than float32
  resolution is one float32 value and is published the same way, in place of
  the previous rejection. The `cellucid` Python package implements the
  identical case.
- Documents the two scopes the gene identifier rules have here, matching the
  `cellucid` Python package. Being drawable is a property of a name the viewer
  shows, so it covers the genes `gene_identifiers` selects; a `var` row left out
  reaches no manifest and is not checked. Distinctness spans the whole `var`,
  because `gene_identifiers` addresses `var` rows by identifier and a repeated
  one names no single row.
- Records whatever `var_gene_id_column` selects faithfully. `cellucid_prepare()`
  performs no symbol lookup and ships no mapping, so the caller decides what a
  gene is called.
- Reports an identifier defect in the same words the `cellucid` Python package
  uses. Each axis passes one singular noun -- `Gene`, `Observation field`,
  `Vector field` -- so the sentence composes as
  `Gene identifier at position 1 ...` instead of
  `Gene identifiers identifier at position 1 ...`, and a repeated key is
  reported as `Gene key 'ACTB' is duplicated.` on both writers. The checks that
  speak about a whole axis add their own plural, so no caller has to guess
  which number a message will need.
- Prints a set of values in a message as a list the caller can paste back into
  the failing call, `c(0, 1, 1)` and `c("score", "n_counts")`. Without the
  boundary, `columns not in obs: a, b. Available columns: x` cannot be read
  back as two lists; with Python's `['score', 'n_counts']` the boundary was
  there but the syntax inside it belonged to the other language. Every message
  that shows a set of values now renders it through one function, so
  `gene_identifiers`, `var_gene_id_column`, and the payload-manifest checks no
  longer each print a set their own way. The `cellucid` Python package prints
  the same sets as Python lists, for the same reason.
- Requires `obs_categorical_dtype`, `dataset_name`, and `dataset_id` as
  arguments with no default, and reports one that was never supplied as
  missing. `cellucid_prepare(dataset_name = "Atlas")` answered
  `dataset_id must be exactly one string.`, which describes a value the caller
  never passed and names nothing valid; it now answers
  `cellucid_prepare() is missing 1 required argument:` followed by
  `dataset_id` and what a valid identifier is. Every argument left out is named
  in that one message, in signature order, as the `cellucid` Python package
  names them in one `TypeError`. A value that *was* supplied and is wrong still
  reports the value: `dataset_id = NULL` remains
  `dataset_id must be exactly one string.`
- Defaults `dataset_description` to `NULL` and accepts `NULL` and `""` alike,
  both publishing `""`. It previously defaulted to `""` and rejected `NULL`,
  while every other optional identity argument in the same signature --
  `source_name`, `source_url`, `source_citation`, `created_at` -- took `NULL`
  for "not supplied". The `cellucid` Python package accepts `None` and `""` for
  it and publishes `""` for both.
- Reports an empty or padded `dataset_id` against the rule it actually breaks.
  `dataset_id = ""` answered
  `dataset_id must be one non-empty string without leading or trailing whitespace.`,
  naming whitespace a value with no characters cannot carry; it now answers
  `dataset_id '' is not a portable identifier. Use 1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a letter or number and not ending with '.'.`
  Nothing that was accepted before is rejected now, and nothing that was
  rejected before is accepted.
- Writes one fixed gzip header for every compressed payload, so a `.gz` no
  longer records the machine that produced it. `gzfile()` leaves the header to
  zlib, which stamps its own build platform into the member's `OS` field --
  `0x03` on Unix, `0x0b` on Windows -- so the same input produced different
  export bytes depending on the operating system the export ran on. Each member
  now carries the ten bytes the `cellucid` Python package writes: no filename,
  a Unix-epoch timestamp, the extra-flags value RFC 1952 assigns to the
  requested level, and `OS = 0xff`, the code that names no operating system.
  The compressed payload, the level the caller chose, and the CRC32 and length
  trailer are untouched, and the deflate stream is still whatever the zlib
  underneath R produces at that level -- so a member matches the Python
  package's byte for byte when the two zlib builds agree, and differs only in
  the deflate bytes when they do not.
- Adds `tests/testthat/test-gzip-header.R`, which reads the ten header bytes of
  a written `.gz` directly at every compression level, and proves the member
  still decompresses to exactly the payload that went in.
- Adds `tests/testthat/test-writer-parity.R`, which holds this writer to the
  half of the two-writer contract it owns: a `vectors/` payload is the scaled
  value rounded to float32 once and not the twice-rounded value, a
  `points_<dim>d.bin` payload is the normalized value rounded once, a
  categorical centroid is measured from coordinates that were never rounded,
  and a `dataset_identity.json` vector field entry carries its keys in the
  documented order.
- Refuses to select row numbers as gene identifiers. A `data.frame` always has
  row names, so the old `is.null(rownames(var))` test could not fire, and
  `var <- data.frame(symbol = c("CD8A", "MS4A1"))` with the default
  `var_gene_id_column = NULL` exported two genes named `"1"` and `"2"` -- the
  automatic row-name sequence that `rownames()` materializes. Those strings are
  unique and drawable, so every later check passed and nothing reported the
  substitution; a wet-lab user then searched `CD8A` in the viewer and found
  nothing. `cellucid_prepare()` now distinguishes automatic row names from row
  names a caller set, and answers `var has only automatic row names, so
  rownames(var) would name the genes '1' to '2'. Set rownames(var) to the gene
  identifiers, or pass var_gene_id_column.` Explicit row names are untouched,
  including ones that happen to read as numbers. The `cellucid` Python package
  already refused the same input, because a default pandas `RangeIndex` yields
  integers where its gene identifiers must be strings.
- Refuses a `data.frame` whose columns carry a class wherever it accepts one in
  place of a matrix. `as.matrix()` keeps the numbers and drops the attribute
  that says what they mean, so a frame of one `units` column in metres and one
  in kilometres became a matrix of equal numbers. `latent_space`,
  `gene_expression`, and the embeddings already rejected such a frame;
  `vector_fields` accepted it and published a 45-degree arrow for data pointing
  0.06 degrees off axis, with no error and nothing downstream able to notice.
  All four axes now apply the one rule. `is.numeric()` answers `FALSE` for
  `Date`, `POSIXct`, `difftime`, and `factor`, so those were already refused;
  what this reaches is every classed numeric with no `is.numeric()` method,
  such as `units::units` and `bit64::integer64`.
- Refuses a classed value for every argument that names an identifier or a
  numeric setting, as the string and logical arguments already did. A classed
  character vector passed as `obs_keys` or `gene_identifiers`, or held in the
  column `var_gene_id_column` selects, reached the check that compares the
  written var manifest against the staged gene names -- `identical()` compares
  attributes -- and failed it with an internal message naming a staging path.
  A classed number passed as `compression`, `var_quantization`,
  `obs_continuous_quantization`, or `centroid_outlier_quantile` reached the
  manifest writer and failed with jsonlite's `No method asJSON S3 class:`,
  which names no argument. Both now report the argument and its rule before any
  file is staged, and the identifier message is
  `<argument> must be a native character vector.`

# cellucid 0.9.0

- Initial GitHub package release.
