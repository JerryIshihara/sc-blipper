# Create conda env

Create and activate conda environment
```bash
conda env create -f environment.yml
conda activate sc-blipper
```

> Note: If you want to give the env a different name: `conda env create -f environment.yml -n my_env`

Note down the install path for the conda env, we will need it later for configuring the pipeline
```bash
echo $CONDA_PREFIX  # On Unix/Linux/macOS
```

Exit conda env
```bash
conda deactivate
```

# Build and run with Apptainer

The unified `environment.yml` can be built once as an OCI image, published to
GHCR, and converted to a shared SIF on the farm. See
[`containers/sc-blipper/README.md`](containers/sc-blipper/README.md) for the
build, registry authentication, upload, pull, and smoke-test commands.

Set the image in the run config, then combine the LSF and Apptainer profiles:

```groovy
params.rn_container = 'file:///lustre/<team>/containers/sc-blipper-poc.sif'
params.preprocess.scvi.container = params.rn_container
params.cnmf_gpu.container = params.rn_container
```

```bash
export NXF_APPTAINER_CACHEDIR=/lustre/<team>/apptainer-cache
sc-blipper cnmf -a -c run.config
```


# Verify GPU dependencies

The shared `environment.yml` includes the CPU, CUDA cNMF, and scVI dependencies.
If you plan to run GPU cNMF or scVI, test the environment on a GPU node:

```python
import jax
import scvi
import torch

print(torch.cuda.is_available())  # Should return True
print(jax.devices())              # Should include a GPU device
print(scvi.__version__)
```
