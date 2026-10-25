"""Read a store cellucid_prepare() wrote with Python's anndata and compare it
with what R put in, value for value.

    python read_with_anndata.py STORE EXPECTED_JSON

Any warning anndata raises while reading is an error. Prints "ok" and exits 0
when every value matches.
"""

import json
import sys
import warnings

warnings.simplefilter("error")

import anndata  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def read(path):
    if path.endswith(".zip"):
        import zarr

        return anndata.read_zarr(zarr.storage.ZipStore(path, mode="r"))
    return anndata.read_zarr(path)


def plain(value):
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [plain(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and value != value:
        # NaN, which the expected JSON writes as null.
        return None
    return value


def main(store, expected_path):
    with open(expected_path, encoding="utf-8") as handle:
        expected = json.load(handle)
    adata = read(store)
    assert list(adata.obs_names) == expected["obs_names"], "obs_names"
    assert list(adata.var_names) == expected["var_names"], "var_names"
    for column, spec in expected["obs"].items():
        values = adata.obs[column]
        if spec["kind"] == "categorical":
            assert isinstance(values.dtype, pd.CategoricalDtype), column
            assert list(values.cat.categories) == spec["categories"], column
            assert values.cat.codes.tolist() == spec["codes"], column
            assert values.cat.ordered == spec["ordered"], column
        else:
            got = [None if pd.isna(item) else plain(item) for item in values]
            assert got == spec["values"], (column, got, spec["values"])
    dense = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    assert adata.X.dtype == np.dtype(expected["X_dtype"]), ("X dtype", adata.X.dtype)
    assert dense.tolist() == expected["X"], "X values"
    if expected["X_format"] == "csc":
        assert adata.X.format == "csc", "X format"
        assert adata.X.nnz == expected["X_nnz"], ("X stored values", adata.X.nnz)
    for key, values in expected["obsm"].items():
        assert adata.obsm[key].dtype == np.float32, key
        assert adata.obsm[key].tolist() == values, key
    graph = adata.obsp["connectivities"]
    assert graph.format == "csr" and graph.dtype == np.float32, "graph"
    assert graph.toarray().tolist() == expected["connectivities"], "graph values"
    assert plain(dict(adata.uns)) == expected["uns"], (plain(dict(adata.uns)), expected["uns"])
    # Reading is over: what the interpreter's own shutdown warns about (some
    # numcodecs releases warn when they clean up) is not about the store.
    warnings.simplefilter("ignore")
    print("ok")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
