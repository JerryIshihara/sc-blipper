FROM condaforge/miniforge3:24.11.3-0@sha256:bcfc04b3562faa8dfd48f3767e2bd7b73a8529f9611b508d42e5bd227db1342d

ARG VCS_REF=unknown
ARG SOURCE_URL=https://github.com/TrynkaLab/sc-blipper

LABEL org.opencontainers.image.title="sc-blipper" \
      org.opencontainers.image.description="Unified CPU, CUDA cNMF, and scVI runtime for sc-blipper" \
      org.opencontainers.image.source="${SOURCE_URL}" \
      org.opencontainers.image.revision="${VCS_REF}"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
        procps \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/sc-blipper-environment.yml

RUN /opt/conda/bin/conda env create \
        --prefix /opt/sc-blipper-env \
        --file /tmp/sc-blipper-environment.yml \
    && /opt/conda/bin/conda clean --all --yes \
    && rm -rf /root/.cache/pip \
    && rm -f /tmp/sc-blipper-environment.yml

COPY --chmod=0755 bin/ /opt/sc-blipper/bin/

ENV PATH="/opt/sc-blipper-env/bin:/opt/sc-blipper/bin:${PATH}" \
    PYTHONNOUSERSITE=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN command -v cnmf \
    && command -v run_preprocess.py \
    && command -v Rscript \
    && python -c "import cnmf, jax, scanpy, scvi, torch; print(torch.__version__)"

WORKDIR /work

CMD ["/bin/bash"]
