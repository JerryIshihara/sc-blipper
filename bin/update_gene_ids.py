#!/usr/bin/env python3
"""
Update gene names in a .h5ad file (optionally using a mapping TSV),
optionally filter genes, and save the result.

Usage:
    python update_gene_ids.py \
        --input input.h5ad \
        --output output.h5ad \
        [--mapping gene_mapping.tsv] \
        [--unmapped unmapped.txt] \
        [--filter-genes keep_genes.txt] \
        [--filter-mode include|exclude]

The mapping TSV should have no header and contain:
    old_id    new_id
"""

import argparse
import os
import pandas as pd
import anndata as ad


def read_h5ad_streaming(input_file, chunk_size=None):
    """Read matrices lazily, reserving 2 GiB or 20% of available host RAM."""
    import h5py
    import math
    import psutil

    available = psutil.virtual_memory().available
    buffer = min(available // 2, max(2 * 1024**3, available // 5))
    budget = available - buffer
    print(f"Streaming RAM budget: {budget / 1024**3:.1f} GiB; buffer: {buffer / 1024**3:.1f} GiB")

    def read_chunked(read, key, elem, *, iospec):
        shape = elem.shape if isinstance(elem, h5py.Dataset) else elem.attrs.get("shape", ())
        if iospec.encoding_type in ("array", "csr_matrix", "csc_matrix") and len(shape) > 1 and all(shape):
            axis = int(iospec.encoding_type == "csc_matrix")
            dtype = elem.dtype if isinstance(elem, h5py.Dataset) else elem["data"].dtype
            # Allow for sparse indices and temporary read/filter/write copies.
            bytes_per_row = 8 * max(dtype.itemsize, 8) * (math.prod(shape) // shape[axis])
            rows = chunk_size or max(1, budget // bytes_per_row)
            chunks = tuple(rows if i == axis else n for i, n in enumerate(shape))
            return ad.experimental.read_elem_as_dask(elem, chunks=chunks)
        return read(elem)

    with h5py.File(input_file, "r") as source:
        return ad.experimental.read_dispatched(source, read_chunked)


def write_h5ad_streaming(adata, output_file):
    """Write one lazy chunk at a time; copy detaches the view without loading X."""
    import dask

    with dask.config.set(scheduler="synchronous"):
        adata.copy().write_h5ad(output_file)


def update_gene_ids(adata, mapping_file, unmapped_file=None):
    """Update gene IDs (var_names) in AnnData using a mapping file."""
    print(f"Reading mapping from {mapping_file} ...")
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["old", "new"])
    mapping_dict = pd.Series(mapping_df["new"].values,
                             index=mapping_df["old"]).to_dict()

    updated_var_names = []
    unmapped_genes = []

    for gene in adata.var_names:
        if gene in mapping_dict:
            updated_var_names.append(mapping_dict[gene])
        else:
            updated_var_names.append(gene)
            unmapped_genes.append(gene)

    adata.var_names = updated_var_names

    if unmapped_genes:
        print(f"Warning: {len(unmapped_genes)} gene IDs could not be mapped and were kept as original.")
        if unmapped_file:
            pd.Series(unmapped_genes).to_csv(unmapped_file, index=False, header=False)
            print(f"Unmapped gene IDs written to: {unmapped_file}")

    return adata, mapping_dict


def filter_genes(adata, filter_file, mapping_dict=None, mode="include"):
    """Filter AnnData object by genes; mode = 'include' or 'exclude'."""
    print(f"Filtering genes using {filter_file} (mode={mode}) ...")
    # Read the filter file and convert to a list
    filter_genes_df = pd.read_csv(filter_file, header=None)
    filter_genes = filter_genes_df[0].astype(str).tolist()  # Access the first column
    print(f"Gene list contains {len(filter_genes)} genes")
    print(f"adata has {adata.shape} cells x genes before filter")
    
    # If mapping provided, update gene IDs in filter list too
    if mapping_dict is not None:
        filter_genes = [mapping_dict.get(g, g) for g in filter_genes]

    if mode == "include":
        keep_mask = adata.var_names.isin(filter_genes)
    elif mode == "exclude":
        keep_mask = ~adata.var_names.isin(filter_genes)
    else:
        raise ValueError("filter-mode must be 'include' or 'exclude'")

    print(f"Keeping {keep_mask.sum()} genes out of {adata.n_vars}")
    adata = adata[:, keep_mask]
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Update and/or filter gene IDs in an AnnData (.h5ad) file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("-o", "--output", required=True, help="Output AnnData file (.h5ad)")
    parser.add_argument("-m", "--mapping", default=None, help="Optional TSV with old_id and new_id (no header)")
    parser.add_argument("-u", "--unmapped", default=None, help="Optional file to write unmapped gene IDs (only if mapping used)")
    parser.add_argument("-f", "--filter-genes", default=None, help="Optional file with list of gene IDs to filter")
    parser.add_argument("--filter-mode", choices=["include", "exclude"], default="include",
                        help="Filter mode: include only listed genes, or exclude them [default: include]")
    parser.add_argument("--streaming", action="store_true", help="Stream matrices from disk (requires anndata[dask]>=0.11)")
    parser.add_argument("--chunk-size", type=int, help="Cells per chunk (CSC: genes); default uses available RAM minus a buffer")

    args = parser.parse_args()
    if args.chunk_size is not None and args.chunk_size < 1:
        parser.error("--chunk-size must be positive")

    # Load AnnData
    print(f"Loading {args.input} ...")
    if args.streaming:
        if os.path.exists(args.output) and os.path.samefile(args.input, args.output):
            parser.error("Streaming input and output must be different files")
        adata = read_h5ad_streaming(args.input, args.chunk_size)
    else:
        adata = ad.read_h5ad(args.input)

    mapping_dict = None
    if args.mapping:
        adata, mapping_dict = update_gene_ids(adata, args.mapping, args.unmapped)
    else:
        print("No mapping file provided. Gene IDs will be kept unchanged.")

    # Optional filtering
    if args.filter_genes:
        adata = filter_genes(adata, args.filter_genes, mapping_dict, mode=args.filter_mode)

    # Save output
    print(f"Saving AnnData to {args.output} ...")
    if args.streaming:
        write_h5ad_streaming(adata, args.output)
    else:
        adata.write_h5ad(args.output)
    print("Done!")


if __name__ == "__main__":
    main()
