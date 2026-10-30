# Sampling workflow smoke test

Requires a cNMF installation containing the current sampling API,
plus the pipeline's normal dependencies.
This is a local CPU wiring test, not a GPU parity or large-data benchmark.

From the repository root, with `cnmf` and `python` from that environment on PATH:

```bash
mkdir -p /tmp/sc-blipper-sampling-test
python tests/nextflow/scripts/make_tiny_counts.py /tmp/sc-blipper-sampling-test/counts.h5ad
touch /tmp/sc-blipper-sampling-test/{NO_TPM,NO_HVG,NO_H5AD}
```

Then run:

```bash
NXF_VER=25.04.6 NXF_SYNTAX_PARSER=v1 nextflow \
  -C tests/nextflow/config/params/sampling_smoke.config \
  run tests/nextflow/sampling_smoke.nf -lib lib \
  --test_counts /tmp/sc-blipper-sampling-test/counts.h5ad \
  --test_tpm /tmp/sc-blipper-sampling-test/NO_TPM \
  --test_hvg /tmp/sc-blipper-sampling-test/NO_HVG \
  --test_no_h5ad /tmp/sc-blipper-sampling-test/NO_H5AD
```

The test uses production processes: one half-size random sample, sampled
prepare, two workers, normal combine, consensus and K-selection at K=2,3.
Four replicates produce 8 and 12 spectra respectively. Both consensus usages
and H5AD exports must have the sampled cell count (100), not the original 200.
There is no full-data prepare, special K-selection gate or streaming override.

To check the original path, add a config override setting
`params.cnmf.sampling.enabled=false`; both Ks should use all 200 cells.

The adapter tests additionally verify sparse count/TPM exports use identical
cell selections, preserve matrix values, recompute sample-only HVGs, and reject
misaligned TPM inputs:

```bash
python -m pytest tests/test_sample_cnmf.py -q
```

## Scope

- Sampling occurs after preprocessing/scVI; it does not reduce scVI memory use.
- HVGs are recomputed from sampled cells through cNMF's chunked HVG selection,
  excluding zero-variance genes. Inherited HVG lists/masks are not reused.
  Only sampled expression goes through prepare and consensus.
- Only `sample_times=1` is supported. Multi-sample consensus and full-cell
  projection are deferred. Workers are not GPU concurrency limits.
- Shared HVGs do not imply shared variance scaling: sampled prepare retains
  cNMF's existing per-sample normalization.
- Sampling defaults off. The runtime/container must contain cNMF sampling
  support; old published images do not gain it from Nextflow configuration.
