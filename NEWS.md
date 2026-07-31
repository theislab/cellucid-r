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
- Prints a set of values in a message as a list, `[0, 1, 1]` and
  `['score', 'n_counts']`, matching the Python package. Without the boundary,
  `columns not in obs: a, b. Available columns: x` cannot be read back as two
  lists.
- Adds `tests/testthat/test-writer-parity.R`, which holds this writer to the
  half of the two-writer contract it owns: a `vectors/` payload is the scaled
  value rounded to float32 once and not the twice-rounded value, a
  `points_<dim>d.bin` payload is the normalized value rounded once, a
  categorical centroid is measured from coordinates that were never rounded,
  and a `dataset_identity.json` vector field entry carries its keys in the
  documented order.
- Adds current multi-platform package checks and a Bootstrap 5 documentation
  site contract.

# cellucid 0.9.0

- Initial GitHub package release.
