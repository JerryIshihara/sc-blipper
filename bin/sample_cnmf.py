#!/usr/bin/env python3
"""Sample once, align TPM, and select HVGs on the sampled cells only."""
import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from cnmf.sampling import sample, _sample_export, _sample_hvgs


def _align_sparse_indices(path):
    """Keep exported CSR/CSC index dtypes consistent for sklearn, without loading X."""
    with h5py.File(path, "r+") as handle:
        groups = []
        handle.visititems(lambda name, obj: groups.append(name)
                          if isinstance(obj, h5py.Group) and
                          obj.attrs.get("encoding-type") in ("csr_matrix", "csc_matrix") else None)
        for name in groups:
            group = handle[name]
            indices, indptr = group["indices"], group["indptr"]
            if indices.dtype == indptr.dtype:
                continue
            dtype = np.int32 if len(indices) <= np.iinfo(np.int32).max and max(group.attrs["shape"]) <= np.iinfo(np.int32).max else np.int64
            for key in ("indices", "indptr"):
                source = group[key]
                if source.dtype == np.dtype(dtype):
                    continue
                target = group.create_dataset(key + "_aligned", shape=source.shape, dtype=dtype)
                for start in range(0, len(source), 1000000):
                    target[start:start + 1000000] = source[start:start + 1000000]
                del group[key]
                group.move(key + "_aligned", key)


def run(args):
    counts = ad.read_h5ad(args.counts, backed="r")
    tpm = None
    try:
        if Path(args.tpm).name != "NO_TPM":
            tpm = ad.read_h5ad(args.tpm, backed="r")
            if not counts.obs_names.equals(tpm.obs_names):
                raise ValueError("Counts and TPM cell IDs/order differ")
            # Full-gene TPM is allowed alongside HVG-only corrected counts.
            if not counts.var_names.isin(tpm.var_names).all():
                raise ValueError("TPM is missing count-matrix genes")
        manifest_path = Path(sample(
            args.counts, args.output, sample_method=args.sample_method,
            sample_fraction=args.sample_fraction, sample_times=args.sample_times,
            sample_seed=args.sample_seed, sample_prefix=args.sample_prefix,
            sample_embedding_key=args.sample_embedding_key,
            sample_numgenes=0,
            sample_export_chunk_size=args.sample_export_chunk_size))
        output = Path(args.output)
        manifest = json.loads(manifest_path.read_text())
        manifest["sampling_stage"] = "after_preprocessing"
        for draw in manifest["samples"]:
            _align_sparse_indices(output / draw["counts"])
            if tpm is None:
                draw["tpm"] = "NO_TPM"
                (output / "NO_TPM").touch()
            else:
                draw["tpm"] = draw["sample_id"] + ".tpm.h5ad"
                indices = np.load(output / draw["cell_indices"])
                _sample_export(tpm, indices, output / draw["tpm"], args.sample_export_chunk_size)
                _align_sparse_indices(output / draw["tpm"])
            # A single sample is an independent cNMF dataset. Do not inherit
            # full-input HVGs or var.highly_variable: genes may be constant here.
            sampled = ad.read_h5ad(output / draw["counts"], backed="r")
            sampled_tpm = None
            try:
                if tpm is not None:
                    sampled_tpm = ad.read_h5ad(output / draw["tpm"], backed="r")
                manifest["hvg"] = _sample_hvgs(
                    sampled, output, draw["sample_id"], args.sample_numgenes,
                    sample_hvg_mode="recompute", tpm=sampled_tpm)
                manifest["shared_hvg"] = manifest["hvg"]["genes_file"]
            finally:
                sampled.file.close()
                if sampled_tpm is not None:
                    sampled_tpm.file.close()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    finally:
        counts.file.close()
        if tpm is not None:
            tpm.file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("counts", "tpm", "hvg", "output"):
        parser.add_argument("--" + key, required=True)
    parser.add_argument("--sample-method", choices=["random", "scsampler"], default="random")
    parser.add_argument("--sample-fraction", type=float, required=True)
    parser.add_argument("--sample-times", type=int, choices=[1], default=1)
    parser.add_argument("--sample-seed", type=int, default=42)
    parser.add_argument("--sample-prefix")
    parser.add_argument("--sample-embedding-key", default="X_pca")
    parser.add_argument("--sample-numgenes", type=int, default=2000)
    parser.add_argument("--sample-export-chunk-size", type=int, default=None)
    run(parser.parse_args())
