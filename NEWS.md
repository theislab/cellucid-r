# cellucid 0.9.1 <!-- CELLUCID_VERSION -->

Version 0.9.1 is the CRAN submission release.

- Exposes `cellucid_prepare()` as the single public preparation entry point.
- Publishes complete prepared generations transactionally and rejects mixed or
  stale output directories.
- Requires exact embedding, observation, gene, vector, and categorical-storage
  inputs before writing any scientific artifact.
- Preserves symmetric positive connectivity weights, exact empty graphs, and
  portable edge-index widths without reinterpreting the source graph.
- Rejects non-finite or non-representable float32 values, duplicate scientific
  identifiers, and non-portable artifact names before publication.
- Adds current multi-platform package checks and a Bootstrap 5 documentation
  site contract.

# cellucid 0.9.0

- Initial GitHub package release.
