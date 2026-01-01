# cellucid 0.99.0 <!-- CELLUCID_VERSION -->

- Initial Bioconductor submission version.
- Adds `cellucid_prepare()` / `prepare()` to export embeddings, metadata, gene expression, and connectivity to the Cellucid viewer format.
- Rejects duplicate gene IDs and filename-sanitization collisions (obs keys / gene IDs) to prevent silent overwrites in exports.
