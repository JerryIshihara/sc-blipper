# sc-blipper container POC

The unified `environment.yml` is the dependency source for both build paths in
this directory. The preferred farm workflow publishes an OCI image to GHCR and
converts it once to a shared, read-only SIF with Apptainer.

## Publish the OCI image to GHCR

Push this branch to GitHub. The branch-specific POC trigger publishes the `poc`,
`tmp-apptainer-poc`, and immutable `sha-<full-commit>` tags automatically:

```bash
git push -u origin tmp/apptainer-poc

gh run list \
  --workflow build-container.yml \
  --branch tmp/apptainer-poc \
  --limit 5

gh run watch <run-id>
```

After this workflow is present on the default branch, it can also be launched
manually with `workflow_dispatch` and an explicit tag. The workflow
authenticates to GHCR with `GITHUB_TOKEN` and builds for `linux/amd64`:

```text
ghcr.io/<github-owner-lowercase>/sc-blipper:poc
ghcr.io/<github-owner-lowercase>/sc-blipper:sha-<full-commit>
```

For a manual push, create a GitHub token with `write:packages`, then run:

```bash
export GHCR_USER=<github-user>
export GHCR_NAMESPACE=<github-owner-lowercase>
export GHCR_TOKEN=<github-token>

printf '%s' "$GHCR_TOKEN" | \
  docker login ghcr.io --username "$GHCR_USER" --password-stdin

docker buildx build \
  --platform linux/amd64 \
  --file containers/sc-blipper/Dockerfile \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  --build-arg SOURCE_URL="$(git remote get-url origin)" \
  --tag "ghcr.io/${GHCR_NAMESPACE}/sc-blipper:poc" \
  --push \
  .
```

New GHCR packages may be private. Either make the package public in its GitHub
package settings or authenticate Apptainer on the farm with a token that has
`read:packages`:

```bash
printf '%s' "$GHCR_TOKEN" | \
  apptainer registry login \
    --username "$GHCR_USER" \
    --password-stdin \
    docker://ghcr.io
```

## Convert the OCI image to SIF on the farm

Use a shared destination and cache that every LSF worker can read:

```bash
export NXF_APPTAINER_CACHEDIR=/lustre/<team>/apptainer-cache
mkdir -p "$NXF_APPTAINER_CACHEDIR" /lustre/<team>/containers

apptainer pull \
  /lustre/<team>/containers/sc-blipper-poc.sif \
  docker://ghcr.io/<github-owner-lowercase>/sc-blipper:poc
```

Prefer the immutable `sha-<full-commit>` tag, or the digest printed by the
workflow, for a production run.

Smoke-test the image:

```bash
IMAGE=/lustre/<team>/containers/sc-blipper-poc.sif

apptainer exec --cleanenv "$IMAGE" cnmf --help

bsub \
  -q gpu \
  -gpu 'num=1:mode=exclusive_process' \
  -o apptainer-gpu-%J.out \
  -e apptainer-gpu-%J.err \
  "apptainer exec --cleanenv --nv $IMAGE python -c \
  'import torch; print(torch.cuda.is_available()); print(torch.cuda.get_device_name())'"
```

## Build and publish a native SIF directly

When fakeroot is available, build from the repository root so `%files` resolves
`environment.yml` and `bin` correctly:

```bash
export APPTAINER_TMPDIR=/local/scratch/$USER/apptainer-tmp
mkdir -p "$APPTAINER_TMPDIR"

apptainer build --fakeroot \
  sc-blipper-poc.sif \
  containers/sc-blipper/sc-blipper.def
```

A native SIF can also be stored in GHCR as an ORAS artifact. The protocol used
to pull must match the protocol used to push:

```bash
printf '%s' "$GHCR_TOKEN" | \
  apptainer registry login \
    --username "$GHCR_USER" \
    --password-stdin \
    oras://ghcr.io

apptainer push --allow-unsigned \
  sc-blipper-poc.sif \
  oras://ghcr.io/${GHCR_NAMESPACE}/sc-blipper-sif:poc

apptainer pull \
  sc-blipper-poc.sif \
  oras://ghcr.io/${GHCR_NAMESPACE}/sc-blipper-sif:poc
```

The OCI route is preferred for this POC because GitHub Actions can build and
publish it without requiring fakeroot on the farm.

This image follows the unified `environment.yml` exactly. It covers the main
sc-blipper, cNMF, GPU cNMF, and scVI processes; the separate CBIIT LDSC runtime
is not included in this POC.
