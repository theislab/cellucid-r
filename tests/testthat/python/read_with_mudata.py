"""Read a MuData store cellucid_prepare() wrote with Python's mudata and print
what it holds as JSON, for test-mudata-interop.R to compare with the R input.

    python read_with_mudata.py STORE

Exits 3 with "skip: <reason>" when mudata is not importable.
"""

import json
import sys

try:
    import mudata
    import numpy as np
except ImportError as error:  # pragma: no cover - depends on the interpreter
    print(f"skip: {error}")
    sys.exit(3)


mdata = mudata.read_zarr(sys.argv[1])


def dense(matrix):
    return np.asarray(matrix.todense() if hasattr(matrix, "todense") else matrix).tolist()


print(json.dumps({
    "n_obs": int(mdata.n_obs),
    "obs_names": list(map(str, mdata.obs_names)),
    "obs_columns": list(map(str, mdata.obs.columns)),
    "mod": list(mdata.mod.keys()),
    "obsm": sorted(map(str, mdata.obsm.keys())),
    "modalities": {
        name: {
            "obs_names": list(map(str, ad.obs_names)),
            "var_names": list(map(str, ad.var_names)),
            "X": dense(ad.X),
            "layers": sorted(key for key in ad.layers.keys() if key is not None),
        }
        for name, ad in mdata.mod.items()
    },
    "obsmap": {name: np.asarray(values).astype(int).tolist() for name, values in mdata.obsmap.items()},
    "cellucid_id": str(mdata.uns["cellucid"]["id"]),
}))
