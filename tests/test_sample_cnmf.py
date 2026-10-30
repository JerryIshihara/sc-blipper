"""Validate the workflow adapter against cNMF's sampling implementation."""
import argparse
import importlib.util
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from cnmf.cnmf import cNMF, compute_tpm, get_highvar_genes_sparse, load_df_from_npz

spec = importlib.util.spec_from_file_location(
    "sample_cnmf", Path(__file__).parents[1] / "bin/sample_cnmf.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


@pytest.mark.parametrize("with_tpm", [False, True])
@pytest.mark.parametrize("with_hvg", [False, True])
def test_aligned_exports(tmp_path, with_tpm, with_hvg):
    x = np.arange(120, dtype=np.float32).reshape(20, 6) + 1
    data = ad.AnnData(sparse.csr_matrix(x))
    data.obs_names = ["cell%d" % i for i in range(20)]
    data.var_names = ["gene%d" % i for i in range(6)]
    counts = tmp_path / "counts.h5ad"
    # Mimic HVG-only corrected counts and full-gene TPM.
    data[:, :3].copy().write_h5ad(counts)
    tpm = tmp_path / ("tpm.h5ad" if with_tpm else "NO_TPM")
    if with_tpm:
        data.write_h5ad(tpm)
    hvg = tmp_path / ("hvg.txt" if with_hvg else "NO_HVG")
    if with_hvg:
        hvg.write_text("\n".join(data.var_names[:3]))
    out = tmp_path / "samples"
    args = argparse.Namespace(counts=str(counts), tpm=str(tpm), hvg=str(hvg),
        output=str(out), sample_method="random", sample_fraction=.25, sample_times=1,
        sample_seed=42, sample_prefix="test", sample_embedding_key="X_pca",
        sample_numgenes=2, sample_export_chunk_size=None)
    adapter.run(args)
    manifest = json.loads((out / "manifest.json").read_text())
    assert len(manifest["samples"]) == 1
    assert manifest["hvg"]["n_input_cells"] == 5
    assert manifest["hvg"]["mode"] == "recompute"
    assert manifest["hvg"]["source"] == "computed"
    assert manifest["shared_hvg"] == manifest["hvg"]["genes_file"]
    assert (out / manifest["shared_hvg"]).is_file()
    for draw in manifest["samples"]:
        indices = np.load(out / draw["cell_indices"])
        subset = ad.read_h5ad(out / draw["counts"])
        assert subset.shape == (5, 3)
        assert sparse.issparse(subset.X)
        assert subset.X.indices.dtype == subset.X.indptr.dtype
        assert subset.obs_names.equals(data.obs_names[indices])
        np.testing.assert_array_equal(subset.X.toarray(), x[indices, :3])
        if with_tpm:
            paired = ad.read_h5ad(out / draw["tpm"])
            assert paired.obs_names.equals(subset.obs_names)
            assert paired.X.indices.dtype == paired.X.indptr.dtype
            np.testing.assert_array_equal(paired.X.toarray(), x[indices])
        else:
            assert draw["tpm"] == "NO_TPM" and (out / "NO_TPM").exists()
        expression = paired.X[:, :3] if with_tpm else compute_tpm(subset).X
        expected, _ = get_highvar_genes_sparse(expression, numgenes=2)
        assert (out / manifest["shared_hvg"]).read_text().splitlines() == subset.var_names[expected.high_var].tolist()


@pytest.mark.parametrize("with_tpm", [False, True])
@pytest.mark.parametrize("streaming", [False, True])
def test_sample_only_hvgs_exclude_zero_variance(tmp_path, with_tpm, streaming):
    # Synthetic full-input HVG gene_0 is expressed only in cells this seeded
    # draw does not select, reproducing the large-data zero-std refit failure.
    values = np.random.default_rng(23).negative_binomial(3, .4, (120, 32)).astype(float) + 1
    seed = int(np.random.SeedSequence(42).generate_state(1)[0])
    selected = np.sort(np.random.default_rng(seed).choice(120, 30, replace=False))
    values[selected, 0] = 0
    data = ad.AnnData(sparse.csr_matrix(values))
    data.var_names = ["gene_%d" % i for i in range(32)]
    data.var["highly_variable"] = [True] + [False] * 31
    counts = tmp_path / "counts.h5ad"
    data.write_h5ad(counts)
    hvg = tmp_path / "full_hvg.txt"
    hvg.write_text("gene_0\n")
    tpm = tmp_path / ("tpm.h5ad" if with_tpm else "NO_TPM")
    if with_tpm:
        compute_tpm(data).write_h5ad(tpm)
    out = tmp_path / "samples"
    args = argparse.Namespace(counts=str(counts), tpm=str(tpm), hvg=str(hvg),
        output=str(out), sample_method="random", sample_fraction=.25, sample_times=1,
        sample_seed=42, sample_prefix="test", sample_embedding_key="X_pca",
        sample_numgenes=40, sample_export_chunk_size=7)
    adapter.run(args)
    manifest = json.loads((out / "manifest.json").read_text())
    draw = manifest["samples"][0]
    np.testing.assert_array_equal(np.load(out / draw["cell_indices"]), selected)
    genes_file = out / manifest["shared_hvg"]
    genes = genes_file.read_text().splitlines()
    assert genes == data.var_names[1:].tolist()
    assert manifest["hvg"]["n_input_cells"] == 30
    run = cNMF(str(tmp_path / "prepared"), "sample")
    run.prepare(str(out / draw["counts"]), [2], n_iter=2, seed=42,
                genes_file=str(genes_file),
                tpm_fn=str(out / draw["tpm"]) if with_tpm else None,
                streaming=streaming)
    normalized = ad.read_h5ad(run.paths["normalized_counts"])
    assert normalized.shape == (30, 31)
    assert np.isfinite(normalized.X.data).all()
    stats = load_df_from_npz(run.paths["tpm_stats"]).loc[genes]
    assert (stats["__std"] > 0).all()


def test_reject_misaligned_tpm(tmp_path):
    data = ad.AnnData(np.ones((6, 3)))
    data.write_h5ad(tmp_path / "counts.h5ad")
    data[::-1].copy().write_h5ad(tmp_path / "tpm.h5ad")
    args = argparse.Namespace(counts=str(tmp_path / "counts.h5ad"), tpm=str(tmp_path / "tpm.h5ad"))
    with pytest.raises(ValueError, match="cell IDs/order"):
        adapter.run(args)
