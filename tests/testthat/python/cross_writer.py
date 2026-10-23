"""Write the input R gave cellucid_prepare() through Python's cellucid.prepare(),
in the Zarr format R wrote, and compare the two stores: the same elements,
encodings, dtypes, shapes, chunks, codecs and fill values, and the same values
in every array (floating values byte for byte, NaN included).

    python cross_writer.py R_STORE INPUT_JSON [R_STORE INPUT_JSON ...]

INPUT_JSON describes the R input in plain JSON (see cross_writer_input() in
tests/testthat/test-cross-writer.R): every dense matrix by rows, every data-frame
column with its R type, every sparse matrix as its compressed-column arrays.
Prints "ok" and exits 0 when every pair of stores agrees; otherwise names the
first difference. Exits 3 with "skip: <reason>" when cellucid, anndata or zarr is
not importable.
"""

import json
import sys
import tempfile
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import zarr
    from scipy import sparse

    import cellucid
except ImportError as error:  # pragma: no cover - depends on the interpreter
    print(f"skip: {error}")
    sys.exit(3)


def matrix(spec):
    if spec["kind"] == "sparse":
        # The compressed-column arrays exactly as R holds them, stored zeros
        # included.
        rows, columns = spec["shape"]
        return sparse.csc_matrix(
            (np.asarray(spec["x"], dtype=np.float64), spec["i"], spec["p"]), shape=(rows, columns)
        )
    dtype = {"double": np.float64, "integer": np.int32, "logical": np.bool_}[spec["type"]]
    return np.asarray(spec["values"], dtype=dtype).reshape(spec["shape"])


def column(spec):
    """A pandas column holding what the R column holds, in R's own types."""
    values = spec["values"]
    kind = spec["type"]
    if kind == "factor":
        return pd.Categorical(values, categories=spec["levels"], ordered=spec["ordered"])
    if kind == "character":
        return np.asarray(values, dtype=object)
    if kind == "double":
        return np.asarray([np.nan if value is None else value for value in values], dtype=np.float64)
    if kind == "integer":
        if None in values:
            return pd.array(values, dtype="Int32")
        return np.asarray(values, dtype=np.int32)
    if kind == "logical":
        if None in values:
            return pd.array(values, dtype="boolean")
        return np.asarray(values, dtype=np.bool_)
    raise ValueError(f"unknown column type {kind}")


def frame(spec, n_rows):
    data = {item["name"]: column(item) for item in spec["columns"]}
    index = spec["index"] if spec["index"] is not None else pd.RangeIndex(n_rows)
    return pd.DataFrame(data, index=index)


def python_arguments(spec):
    n_cells = next(iter(spec["embeddings"].values()))["shape"][0]
    arguments = dict(spec["arguments"])
    arguments["embeddings"] = {key: matrix(value) for key, value in spec["embeddings"].items()}
    arguments["obs"] = frame(spec["obs"], n_cells)
    if spec.get("gene_expression") is not None:
        expression = matrix(spec["gene_expression"])
        arguments["gene_expression"] = expression
        arguments["var"] = frame(spec["var"], expression.shape[1])
    if spec.get("connectivities") is not None:
        arguments["connectivities"] = matrix(spec["connectivities"])
    if spec.get("vector_fields") is not None:
        arguments["vector_fields"] = {
            field: {key: matrix(value) for key, value in bindings.items()}
            for field, bindings in spec["vector_fields"].items()
        }
    if spec.get("embedding_defaults") is not None:
        arguments["embedding_defaults"] = {
            slot: value if isinstance(value, str) else (value["key"], value["columns"])
            for slot, value in spec["embedding_defaults"].items()
        }
    return arguments


def zarr_format(store):
    return 3 if Path(store, "zarr.json").is_file() else 2


def metadata(store):
    """Every metadata document of the store, as on disk, by key.

    Zarr format 2: every document ``.zmetadata`` lists, read from its own
    file (zarr-python's consolidated copy of a nested group carries an extra,
    empty ``consolidated_metadata`` member that no reader of a format 2 store
    uses). Zarr format 3: every node's ``zarr.json``, and the root's
    consolidated metadata as one more document.
    """
    if zarr_format(store) == 2:
        with open(Path(store, ".zmetadata"), encoding="utf-8") as handle:
            listed = json.load(handle)["metadata"]
        documents = {}
        for key in listed:
            with open(Path(store, key), encoding="utf-8") as handle:
                documents[key] = json.load(handle)
        return documents
    with open(Path(store, "zarr.json"), encoding="utf-8") as handle:
        root = json.load(handle)
    consolidated = root.pop("consolidated_metadata")
    documents = {"zarr.json": root, "consolidated_metadata": consolidated}
    for path in consolidated["metadata"]:
        with open(Path(store, path, "zarr.json"), encoding="utf-8") as handle:
            documents[f"{path}/zarr.json"] = json.load(handle)
    return documents


def array_paths(documents):
    """The path of every array the documents describe."""
    paths = []
    for key, document in documents.items():
        if key.endswith(".zarray"):
            paths.append(key[: -len("/.zarray")] if "/" in key else "")
        elif key.endswith("zarr.json") and document.get("node_type") == "array":
            paths.append(key[: -len("/zarr.json")])
    return sorted(paths)


def nullable_values(path, root):
    """Whether ``path`` is the values array of a nullable column: unspecified where the mask is set."""
    if not path.endswith("/values"):
        return False
    parent = path.rsplit("/", 1)[0]
    return str(root[parent].attrs.get("encoding-type", "")).startswith("nullable-")


def fail(message):
    print(message)
    sys.exit(1)


def compare(r_store, python_store):
    if zarr_format(r_store) != zarr_format(python_store):
        fail(f"R wrote Zarr format {zarr_format(r_store)}, Python {zarr_format(python_store)}")
    r_documents = metadata(r_store)
    python_documents = metadata(python_store)
    if set(r_documents) != set(python_documents):
        fail(
            "different elements: only R "
            f"{sorted(set(r_documents) - set(python_documents))}, only Python "
            f"{sorted(set(python_documents) - set(r_documents))}"
        )
    for key in sorted(r_documents):
        if r_documents[key] != python_documents[key]:
            fail(f"{key} differs: R {r_documents[key]} Python {python_documents[key]}")
    r_root = zarr.open_group(r_store, mode="r")
    python_root = zarr.open_group(python_store, mode="r")
    for path in array_paths(r_documents):
        r_values = np.asarray(r_root[path][...] if path else r_root[...])
        python_values = np.asarray(python_root[path][...] if path else python_root[...])
        if nullable_values(path, r_root):
            mask = np.asarray(r_root[path.rsplit("/", 1)[0] + "/mask"][...])
            r_values = r_values[~mask]
            python_values = python_values[~mask]
        if r_values.dtype != python_values.dtype or r_values.shape != python_values.shape:
            fail(f"{path}: R {r_values.dtype}{r_values.shape}, Python {python_values.dtype}{python_values.shape}")
        # Floating values compare by their bytes: NaN equals no value, and
        # which NaN (and which zero) is stored is part of the output.
        same = (
            r_values.tobytes() == python_values.tobytes()
            if r_values.dtype.kind == "f"
            else np.array_equal(r_values, python_values)
        )
        if not same:
            fail(f"{path} values differ: R {r_values.tolist()} Python {python_values.tolist()}")


def main(pairs):
    for r_store, input_path in pairs:
        with open(input_path, encoding="utf-8") as handle:
            spec = json.load(handle)
        with tempfile.TemporaryDirectory() as folder:
            python_store = Path(folder, "python.zarr")
            cellucid.prepare(
                python_store, format="zarr", zarr_format=zarr_format(r_store), **python_arguments(spec)
            )
            compare(r_store, str(python_store))
    print("ok")


if __name__ == "__main__":
    arguments = sys.argv[1:]
    main(list(zip(arguments[::2], arguments[1::2], strict=True)))
