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
- Rejects non-finite or non-representable float32 values, duplicate scientific
  identifiers, and non-portable artifact names before publication.
- Names the rejected identifier axis when an artifact name is not portable or is
  a Windows device name, so gene identifiers, `obs_keys`, observation field
  keys, and vector field ids no longer all report `Filename component`.
- Requires every string the viewer prints verbatim to read on screen as the
  value it stores. A string category label, `dataset_name`,
  `dataset_description`, `source_name`, `source_url`, and `source_citation` are
  rejected when they carry a control character, one of the zero-width
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
  float32 value domain and rejects ranges whose variation would disappear.
- Documents the two scopes the gene identifier rules have always had here, now
  that the `cellucid` Python package matches them. The portable-filename and
  case-insensitive collision rules are about the file a gene is written to, so
  they cover the genes `gene_identifiers` selects; a `var` row left out is
  written to no path and is not checked. The rule that is not about paths spans
  the whole `var`: every gene identifier must be a non-empty string and must be
  distinct, because `gene_identifiers` addresses `var` rows by identifier.
- Adds current multi-platform package checks and a Bootstrap 5 documentation
  site contract.

# cellucid 0.9.0

- Initial GitHub package release.
