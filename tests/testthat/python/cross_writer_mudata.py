"""Write a MuData described in JSON through Python's cellucid.prepare(), and
compare the MuData Zarr store with one another writer made from the same input.

    python cross_writer_mudata.py OTHER_STORE INPUT_JSON [OTHER_STORE INPUT_JSON ...]

It is the Python half of the MuData cross-writer test: cellucid-r's
``cellucid_prepare(data = <multi-assay object>)`` writes ``OTHER_STORE`` (a
MuData Zarr folder) and describes the very same input in ``INPUT_JSON``; this
script builds that MuData, writes it with ``cellucid.prepare(data=mdata)`` in
the Zarr format of ``OTHER_STORE``, and compares the two stores: the same
groups and arrays, attributes, dtypes, shapes, chunks, codecs, fill values,
and every value (floating values byte for byte, NaN included). Prints "ok" and
exits 0 when every pair agrees; otherwise names the first difference. Exits 3
with "skip: <reason>" when cellucid, mudata or zarr is not importable.

INPUT_JSON::

    {
      "arguments": {"dataset_id": ..., "dataset_name": ..., "modality": ...},
      "obs": FRAME,                       # the global cells, in global order
      "obsm": {"X_wnn_umap": MATRIX},     # shared embeddings
      "obsp": {"connectivities": MATRIX}, # the shared neighbour graph, if any
      "modalities": [                     # in mod-order
        {"name": "rna", "X": MATRIX, "obs": FRAME, "var": FRAME, "obsm": {...}}
      ]
    }

A MATRIX is ``{"kind": "sparse", "shape": [rows, columns], "p", "i", "x"}``
(compressed columns, cells x features) or ``{"kind": "dense", "type":
"double"|"integer", "shape", "values"}`` (values by rows). A FRAME is
``{"index": [names], "columns": [{"name", "type", "values", "levels"?,
"ordered"?}]}`` with R's column types (``factor``, ``character``, ``double``,
``integer``, ``logical``), as the AnnData cross-writer input
(``tests/testthat/python/cross_writer.py`` in cellucid-r) writes them.

The same file is kept, byte for byte, in cellucid-r as
``tests/testthat/python/cross_writer_mudata.py``; both repositories' tests
check the two copies are identical.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from typing import Any

try:
    import anndata as ad
    import mudata as md
    import numpy as np
    import pandas as pd
    import zarr
    from scipy import sparse

    import cellucid
except ImportError as error:  # pragma: no cover - depends on the interpreter
    print(f"skip: {error}")
    sys.exit(3)


def matrix(spec: dict[str, Any]) -> Any:
    if spec["kind"] == "sparse":
        rows, columns = spec["shape"]
        return sparse.csc_matrix(
            (np.asarray(spec["x"], dtype=np.float64), spec["i"], spec["p"]), shape=(rows, columns)
        )
    dtype = {"double": np.float64, "integer": np.int32, "logical": np.bool_}[spec["type"]]
    return np.asarray(spec["values"], dtype=dtype).reshape(spec["shape"])


def column(spec: dict[str, Any]) -> Any:
    values = spec["values"]
    kind = spec["type"]
    if kind == "factor":
        return pd.Categorical(values, categories=spec["levels"], ordered=spec["ordered"])
    if kind == "character":
        return np.asarray(values, dtype=object)
    if kind == "double":
        return np.asarray([np.nan if value is None else value for value in values], dtype=np.float64)
    if kind == "integer":
        return pd.array(values, dtype="Int32") if None in values else np.asarray(values, dtype=np.int32)
    if kind == "logical":
        return pd.array(values, dtype="boolean") if None in values else np.asarray(values, dtype=bool)
    raise ValueError(f"unknown column type {kind!r}")


def frame(spec: dict[str, Any]) -> pd.DataFrame:
    data = {item["name"]: column(item) for item in spec["columns"]}
    index = pd.Index([str(name) for name in spec["index"]], dtype=object)
    return pd.DataFrame(data, index=index)


def mudata_from_spec(spec: dict[str, Any]) -> Any:
    """The MuData ``INPUT_JSON`` describes."""
    modalities: dict[str, Any] = {}
    for item in spec["modalities"]:
        adata = ad.AnnData(X=matrix(item["X"]), obs=frame(item["obs"]), var=frame(item["var"]))
        for key, value in item.get("obsm", {}).items():
            adata.obsm[key] = matrix(value)
        modalities[item["name"]] = adata
    mdata = md.MuData(modalities)
    obs = frame(spec["obs"])
    if list(obs.index) != list(mdata.obs_names):
        raise ValueError(
            "The input's global cells are not in mudata's global order "
            f"({list(obs.index)[:5]} against {list(mdata.obs_names)[:5]})"
        )
    for name in obs.columns:
        mdata.obs[name] = obs[name].to_numpy() if not isinstance(obs[name].dtype, pd.CategoricalDtype) else obs[name].values
    for key, value in spec.get("obsm", {}).items():
        mdata.obsm[key] = matrix(value)
    for key, value in spec.get("obsp", {}).items():
        mdata.obsp[key] = sparse.csr_matrix(matrix(value))
    return mdata


def _codecs(array: Any) -> list[Any]:
    metadata = array.metadata.to_dict()
    return [metadata.get(key) for key in ("compressor", "filters", "codecs", "compressors", "dtype", "data_type", "fill_value", "chunks", "chunk_grid")]


def compare_stores(left: Any, right: Any, path: str = "") -> str | None:
    """The first difference between two Zarr groups, or None."""
    if dict(left.attrs) != dict(right.attrs):
        return f"{path or '/'}: attributes {dict(left.attrs)} != {dict(right.attrs)}"
    left_keys, right_keys = sorted(left.keys()), sorted(right.keys())
    if left_keys != right_keys:
        return f"{path or '/'}: members {left_keys} != {right_keys}"
    for key in left_keys:
        a, b = left[key], right[key]
        where = f"{path}/{key}"
        if isinstance(a, zarr.Group) != isinstance(b, zarr.Group):
            return f"{where}: a group in one store, an array in the other"
        if isinstance(a, zarr.Group):
            difference = compare_stores(a, b, where)
            if difference is not None:
                return difference
            continue
        if dict(a.attrs) != dict(b.attrs):
            return f"{where}: attributes {dict(a.attrs)} != {dict(b.attrs)}"
        if _codecs(a) != _codecs(b):
            return f"{where}: layout {_codecs(a)} != {_codecs(b)}"
        values_a, values_b = np.asarray(a[...]), np.asarray(b[...])
        if values_a.dtype != values_b.dtype:
            same = False
        elif values_a.dtype.kind == "f":
            same = values_a.tobytes() == values_b.tobytes()
        else:
            same = np.array_equal(values_a, values_b)
        if not same:
            return f"{where}: values differ"
    return None


def compare(other_store: Path, input_json: Path) -> str | None:
    spec = json.loads(input_json.read_text(encoding="utf-8"))
    other = zarr.open_group(str(other_store), mode="r")
    zarr_format = other.metadata.zarr_format
    with tempfile.TemporaryDirectory() as folder:
        mine = cellucid.prepare(
            Path(folder) / "python.zarr",
            format="zarr",
            zarr_format=zarr_format,
            data=mudata_from_spec(spec),
            **spec.get("arguments", {}),
        )
        return compare_stores(zarr.open_group(str(mine), mode="r"), other)


def main(arguments: list[str]) -> int:
    if not arguments or len(arguments) % 2:
        print(__doc__)
        return 2
    for other, described in zip(arguments[::2], arguments[1::2], strict=True):
        difference = compare(Path(other), Path(described))
        if difference is not None:
            print(f"{other}: {difference}")
            return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
