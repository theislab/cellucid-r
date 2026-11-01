"""Serve, through cellucid-python's served store, the object R served, and
print the SHA-256 of every key of the arrays both serve.

    python served_store.py INPUT_JSON KEY [KEY ...]

INPUT_JSON describes the object (see test-served-cross-language.R): the cell
and gene names, the cells x genes matrix as its compressed-column arrays (or
as coordinates, float64, duplicates and stored zeros kept), the UMAP by
rows, and optionally a symmetric neighbour graph (``graph``: its
compressed-column arrays, float32 weights), stored as
``obsp['connectivities']``. A matrix the store refuses prints ``{"__refused__": message}``. Prints one JSON object: each array's .zarray as its
metadata, and the SHA-256 of every chunk of the array, and under
``__chunk_elements__`` the chunk sizes the Python served store uses
(``cellucid.prepare_data._store``: per cell and whole; a served expression
matrix's entries, ``cellucid._served_store``). Exits 3 with
"skip: <reason>" when cellucid or anndata is not importable.
"""

import hashlib
import json
import sys

try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    from cellucid._served_store import SERVED_ENTRY_CHUNK_ELEMENTS, PriorityStore
    from cellucid.prepare_data._store import CELL_CHUNK_ELEMENTS, WHOLE_CHUNK_ELEMENTS
except ImportError as error:  # pragma: no cover - depends on the interpreter
    print(f"skip: {error}")
    sys.exit(3)


def main(path: str, arrays: list[str]) -> None:
    spec = json.loads(open(path, encoding="utf-8").read())
    n, m = len(spec["cells"]), len(spec["genes"])
    if "rows" in spec:
        # Coordinates as given, float64: duplicates, stored zeros and values
        # float32 rounds kept as they are (a CSR matrix scipy would hold).
        rows = np.asarray(spec["rows"], dtype=np.int64)
        order = np.argsort(rows, kind="stable")
        x = sparse.csr_matrix(
            (
                np.asarray([float(v) if v is not None else np.nan for v in spec["values"]], dtype=np.float64)[order],
                np.asarray(spec["cols"], dtype=np.int32)[order],
                np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=n))]),
            ),
            shape=(n, m),
        )
    else:
        x = sparse.csc_matrix(
            (np.asarray(spec["x"], dtype=np.float32), np.asarray(spec["i"], dtype=np.int32), np.asarray(spec["p"], dtype=np.int32)),
            shape=(n, m),
        ).tocsr()
    obsp = {}
    if "graph" in spec:
        graph = spec["graph"]
        # Symmetric: its compressed-column arrays are its compressed-row ones.
        obsp["connectivities"] = sparse.csr_matrix(
            (np.asarray(graph["x"], dtype=np.float32), np.asarray(graph["i"], dtype=np.int32), np.asarray(graph["p"], dtype=np.int32)),
            shape=(n, n),
        )
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=pd.Index(spec["cells"], dtype=object)),
        var=pd.DataFrame(index=pd.Index(spec["genes"], dtype=object)),
        obsm={"X_umap": np.asarray(spec["umap"], dtype=np.float32).reshape(n, 2)},
        obsp=obsp,
    )
    try:
        store = PriorityStore.of_anndata(adata)
    except ValueError as error:
        # Refused (a value that is not finite), as prepare() refuses it.
        print(json.dumps({"__refused__": str(error)}))
        return
    try:
        # Its expression matrices are laid out by the store's worker: wait for them.
        store.laid_out()
        keys = list(store.stored_keys())
        out = {}
        for array in arrays:
            for key in keys:
                if key == f"{array}/.zarray" or (key.startswith(f"{array}/") and not key.rsplit("/", 1)[1].startswith(".")):
                    body = store.get(key)
                    assert body is not None, key
                    out[key] = json.loads(body) if key.endswith(".zarray") else hashlib.sha256(body).hexdigest()
        out["__chunk_elements__"] = {"cell": CELL_CHUNK_ELEMENTS, "served": SERVED_ENTRY_CHUNK_ELEMENTS, "whole": WHOLE_CHUNK_ELEMENTS}
        print(json.dumps(out))
    finally:
        store.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
