# cellucid 0.9.1 <!-- CELLUCID_VERSION -->

Version 0.9.1 is the CRAN submission release.

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
  scientific identifiers before publication.
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
