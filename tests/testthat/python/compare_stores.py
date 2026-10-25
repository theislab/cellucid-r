"""Read two Zarr stores with zarr-python and compare them array by array.

    python compare_stores.py PREPARED SERVED

PREPARED is a store cellucid_prepare() wrote; SERVED is the store
cellucid_view() serves from memory, written key by key to a folder. Both must
hold the same groups and arrays, with the same attributes, and every array
the same dtype, shape, chunk shape and values (NaN equal to NaN). The served
store declares no compressor and no filter but its text codec; the prepared
one compresses. Prints one line per difference and exits 1 when there is
any, else prints "ok N" (the number of arrays compared) and exits 0.
"""

import sys
import warnings

warnings.simplefilter("error")

import numpy as np  # noqa: E402
import zarr  # noqa: E402


def nodes(root):
    found = {}

    def walk(group, prefix):
        for name, member in group.members():
            path = f"{prefix}{name}"
            found[path] = member
            if isinstance(member, zarr.Group):
                walk(member, f"{path}/")

    walk(root, "")
    return found


def same_values(left, right):
    left, right = np.asarray(left), np.asarray(right)
    if left.dtype != right.dtype:
        return False
    if left.dtype.kind == "f":
        return np.array_equal(left, right, equal_nan=True)
    return np.array_equal(left, right)


def main(prepared_path, served_path):
    prepared = zarr.open_group(prepared_path, mode="r")
    served = zarr.open_group(served_path, mode="r")
    differences = []
    if dict(prepared.attrs) != dict(served.attrs):
        differences.append("/: attributes differ")
    left, right = nodes(prepared), nodes(served)
    if set(left) != set(right):
        differences.append(f"nodes differ: {sorted(set(left) ^ set(right))}")
    arrays = 0
    for path in sorted(set(left) & set(right)):
        a, b = left[path], right[path]
        if type(a) is not type(b):
            differences.append(f"{path}: a {type(a).__name__} and a {type(b).__name__}")
            continue
        if dict(a.attrs) != dict(b.attrs):
            differences.append(f"{path}: attributes differ")
        if isinstance(a, zarr.Group):
            continue
        arrays += 1
        if b.compressors:
            differences.append(f"{path}: the served array declares a compressor {b.compressors}")
        if a.dtype != b.dtype or a.shape != b.shape or a.chunks != b.chunks:
            differences.append(
                f"{path}: {a.dtype} {a.shape} chunks {a.chunks} and {b.dtype} {b.shape} chunks {b.chunks}"
            )
            continue
        if not same_values(a[...], b[...]):
            differences.append(f"{path}: values differ")
    for difference in differences:
        print(difference)
    if differences:
        return 1
    print(f"ok {arrays}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
