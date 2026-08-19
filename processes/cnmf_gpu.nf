#!/usr/bin/env nextflow

// Prepare itself does not use a GPU. This variant runs on the ordinary cNMF
// queue while using the GPU-capable environment and cNMF's GPU-aware prepare
// adapter, which persists the requested solver for later GPU stages.
process cnmf_prepare_gpu {
    label params.cnmf.label
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.prepare?.container ?: params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.prepare?.conda ?: params.cnmf_gpu?.conda ?: params.rn_conda}"

    input:
        tuple val(id), path(file)
        path(tpm)
        path(hvg)
    output:
        tuple val(id), path("${id}")

    script:
        def gpu_args = CnmfGpuConfigurer.process(params.cnmf_gpu)
        def seed = Math.round(Math.random() * 10000).toInteger()
        cmd =
        """
        cnmf prepare \
        --output-dir ./ \
        --name ${id} \
        -c ${file} \
        -k ${params.cnmf.k.split(",").join(' ')} \
        --n-iter ${params.cnmf.n_iter} \
        --numgenes ${params.preprocess.n_variable}\
        """

        if (params.cnmf.seed != null) {
            cmd += " --seed ${params.cnmf.seed}"
        } else {
            cmd += " --seed ${seed}"
        }

        if (tpm.getFileName().toString() != "NO_TPM") {
            cmd += " --tpm ${tpm}"
        }

        if (hvg.getFileName().toString() != "NO_HVG") {
            cmd += " --genes-file ${hvg}"
        }
        
        cmd += " --engine gpu ${gpu_args.prepare}"

        cmd
}

process cnmf_factorize_gpu {
    label "${params.cnmf_gpu?.label ?: 'gpu_medium'}"
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.conda ?: params.rn_conda}"

    input:
        tuple val(id), path(file, name: "tmp/*"), val(worker_index)
        val n_workers
    output:
        //tuple val(id), path("${id}/cnmf_tmp/*.df.npz")
        path("${id}/cnmf_tmp/*.k_*.iter_*.df.npz", emit: files)

    script:
        def gpu_args = CnmfGpuConfigurer.process(params.cnmf_gpu)

        cmd =
        """
        mkdir -p ${id}/cnmf_tmp

        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/

        cnmf factorize \
        --output-dir ./ \
        --name ${id} \
        --worker-index ${worker_index} \
        --total-workers ${n_workers} \
        --engine gpu ${gpu_args.factorize}
        """

        cmd
}

process cnmf_consensus_gpu {
    label "${params.cnmf_gpu?.label ?: 'gpu_medium'}"
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.conda ?: params.rn_conda}"
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_${k}", mode: 'symlink', saveAs: { filepath ->
        def pathStr = filepath.toString()
        def newPathStr = pathStr.replaceFirst("^${id}/", "")
        return newPathStr
    }

    input:
        tuple val(id), path(file, name: "tmp/*"), val(k)
        path(merged, name: "merged/*")
        path(h5ad)
    output:
        tuple val(k), path("${id}/${id}.gene_spectra_score*"), emit: spectra_score
        tuple val(k), path("${id}/${id}.gene_spectra_tpm*"), emit: spectra_tpm
        tuple val(k), path("${id}/${id}.spectra.k_*"), emit: spectra_k
        tuple val(k), path("${id}/${id}.usages.k_*"), emit: usages_k
        tuple val(k), path("${id}/${id}.starcat_spectra.k_*"), emit: starcat_spectra_k
        tuple val(k), path("${id}/${id}.*.png"), emit: plots
        tuple val(k), path("${id}/${id}.*.h5ad"), optional: true, emit: h5ad
        tuple val(k), path("${id}/cnmf_tmp/${id}.local_density_cache.k_${k}.merged.df.npz"), path("${id}/cnmf_tmp/${id}.spectra.k_${k}.merged.df.npz"),  path("${id}/cnmf_tmp/${id}.nmf_idvrun_params.yaml"),  path("${id}/cnmf_tmp/${id}.norm_counts.h5ad"), emit: qc_input
    script:

        String local_dens = params.cnmf.local_density.toString().replace('.', '_')
        def gpu_args = CnmfGpuConfigurer.process(params.cnmf_gpu)

        cmd =
        """
        mkdir -p ${id}/cnmf_tmp

        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/
        ln -s \$(pwd)/merged/* ${id}/cnmf_tmp/

        cnmf consensus \
        --output-dir ./ \
        --name ${id} \
        --components ${k} \
        --local-density-threshold ${params.cnmf.local_density} \
        --show-clustering \
        --engine gpu ${gpu_args.consensus}

        find ${id} -maxdepth 1 -type f -name "*.txt" ! -xtype l -exec gzip {} \\;
        """

        if (h5ad.getFileName().toString() != "NO_H5AD") {
            cmd +=
            """
            cnmf_to_h5ad.py \
            --spectra ${id}/${id}.gene_spectra_score.*.txt.gz \
            --usage ${id}/${id}.usages.k_*.txt.gz \
            --obs ${h5ad} \
            --output ${id}/${id}.k_${k}.dt_${local_dens}.h5ad
            """
        }

        cmd
}
