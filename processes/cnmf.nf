#!/usr/bin/env nextflow

process cnmf_pre_process {
    label params.cnmf.label_high
    scratch params.rn_scratch
    
    container params.rn_container
    containerOptions Configurer.process(params).cnmf.preprocess.containerOptions
    conda params.rn_conda

    publishDir "$params.rn_publish_dir/h5ad/cnmf/${params.rn_runname}", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
    output:
        path("*.Corrected.HVG.Varnorm.h5ad", emit: corrected)
        path("*.TP10K.h5ad", emit: normalized)
        path("*.Corrected.HVGs.txt", emit: hvg)
        tuple val(id), path("*.Corrected.HVG.Varnorm.h5ad"), emit: preprocessed
        
    script:
        cmd = "run_preprocess.py -i ${file} -o ${id}"
        
        if (params.convert.harmony_vars != null) {
            cmd += " --harmony_vars ${params.convert.harmony_vars.split(',').join(' ')}"
        }
        
        if (params.cnmf.seed != null) {
            cmd += " --seed ${params.cnmf.seed}"
        }
        
        if (params.preprocess.n_variable != null) {
            cmd += " --n_variable ${params.preprocess.n_variable}"
        }
        
        if (params.cnmf.feature_type_col != null) {
            cmd += " --feature_type_col ${params.cnmf.feature_type_col}"
        }
        
        cmd
}

process cnmf_prepare {
    label Configurer.process(params).cnmf.prepare.label
    scratch Configurer.process(params).cnmf.prepare.scratch

    // Prepare does not request a GPU, but GPU-enabled runs need the cNMF GPU
    // package so the selected solver and loss can be persisted for later stages.
    container Configurer.process(params).cnmf.prepare.container
    conda Configurer.process(params).cnmf.prepare.conda

    // This doesn't really need to be published
    //publishDir "$params.rn_publish_dir/cnmf/prepare/", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
        path(tpm)
        path(hvg)
    output:
        tuple val(id), path("${id}")
        
    script:
        def config = Configurer.process(params).cnmf
        def seed =  Math.round( Math.random() * 10000 ).toInteger()
        cmd = 
        """
        cnmf prepare \
        --output-dir ./ \
        --name $id \
        -c $file \
        -k ${config.common.k.join(' ')} \
        --n-iter ${config.common.n_iter} \
        --numgenes ${config.preprocess.n_variable}\
        """
        
        if (config.common.seed != null) {
            cmd += " --seed ${config.common.seed}"
        } else {
            cmd += " --seed ${seed}"
        }
        
        // Add the TPM matrix/h5ad to calculate the gene scores
        if (tpm.getFileName().toString() != "NO_TPM") {
            cmd += " --tpm $tpm"
        } 
        
        // Add the variable genes, if providing TPM this must also be provided
        if (hvg.getFileName().toString() != "NO_HVG") {
            cmd += " --genes-file $hvg"
        }

        if (config.prepare.args) {
            cmd += " ${config.prepare.args}"
        }

        cmd
}

process cnmf_factorize {
    label Configurer.process(params).cnmf.factorize.label
    scratch Configurer.process(params).cnmf.factorize.scratch

    container Configurer.process(params).cnmf.factorize.container
    containerOptions Configurer.process(params).cnmf.factorize.containerOptions
    conda Configurer.process(params).cnmf.factorize.conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    //publishDir "$params.rn_publish_dir/cnmf/${id}/factorize", mode: 'symlink'
    
    input:
        tuple val(id), path(file, name: "tmp/*"), val(worker_index)
        val n_workers
    output:
        //tuple val(id), path("${id}/cnmf_tmp/*.df.npz")
        path("${id}/cnmf_tmp/*.k_*.iter_*.df.npz", emit: files)
        
    script:
        def config = Configurer.process(params).cnmf
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir    
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/

        # Now we can run the process, so the output is written to the ${id}/cnmf_tmp in the workdir, 
        # not the symlinked one
        cnmf factorize \
        --output-dir ./ \
        --name ${id} \
        --worker-index ${worker_index} \
        --total-workers ${n_workers}
        """.stripIndent().trim()

        if (config.factorize.args) {
            cmd += " ${config.factorize.args}"
        }
    
        cmd
}

process cnmf_combine {
    label params.cnmf.label_high
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    //publishDir "$params.rn_publish_dir/cnmf/${id}/combine", mode: 'symlink'
    
    input:
        tuple val(id), path(file, name: "tmp/*")
        path(factorize, name: "factorized/*")
    output:
        //tuple val(id), path("${id}")
        path("${id}/cnmf_tmp/*.k_*.merged.df.npz", emit: files)
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/
        ln -s \$(pwd)/factorized/* ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        cnmf combine \
        --output-dir ./ \
        --name ${id}
        """
    
        cmd
}

process cnmf_kselection {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_selection", mode: 'symlink'
    
    input:
        tuple val(id), path(file, name: "tmp/*")
        path(merged, name: "merged/*")
    output:
        path("${id}.k_selection*")
        
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/
        ln -s \$(pwd)/merged/* ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        
        # K selection plot for cNMF
        cnmf k_selection_plot --output-dir ./ --name ${id}
        
        # Symlink for staging
        ln -s ${id}/${id}.k_selection* ./
        """
    
        cmd
}

process cnmf_consensus {
    label params.cnmf.label_high
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_${k}", mode: 'symlink', saveAs: { filepath ->
        // filepath is a Path object, convert to string
        def pathStr = filepath.toString()
        // Remove the first folder (id/)
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
        tuple val(k), path("${id}/${id}*.annotation.tsv"), path("${id}/${id}*.edist.tsv"), emit: qc_metrics
        path("${id}/${id}*.edist.pdf"), emit: qc_plots
    script:
    
        String local_dens = params.cnmf.local_density.toString().replace('.', '_')
    
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/
        ln -s \$(pwd)/merged/* ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        
        # Consensus for cNMF
        cnmf consensus \
        --output-dir ./ \
        --name ${id} \
        --components ${k} \
        --local-density-threshold ${params.cnmf.local_density} \
        --show-clustering
        
        # Zip the textfiles
        find ${id} -maxdepth 1 -type f -name "*.txt" ! -xtype l -exec gzip {} \\;
        """
    
        if (h5ad.getFileName().toString() != "NO_H5AD") {
            cmd +=
            """
            # Convert to h5ad file
            cnmf_to_h5ad.py \
            --spectra ${id}/${id}.gene_spectra_score.*.txt.gz \
            --usage ${id}/${id}.usages.k_*.txt.gz \
            --obs ${h5ad} \
            --output ${id}/${id}.k_${k}.dt_${local_dens}.h5ad
            """
        }
    
    
        cmd
}


process cnmf_ktree {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_selection/", mode: 'symlink'
    
    input:
        path(files)
        val(qry_string)
        val(id)
    output:
        path("${id}*.pdf", emit: plots)
        path("${id}*.tsv", emit: edge_list)
    script:
        cmd = 
        """
        plot_k_tree.r \
        --qry "${qry_string}" \
        --threshold ${params.cnmf.ktree_threshold} \
        --mode ${params.cnmf.ktree_mode} \
        --output $id
        """
    
        cmd
}

process cnmf_summarize {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/${k}/final", mode: 'symlink'
    
    input:
        val(id)
        tuple val(k), path(spectra_file), path(enrich_file), path(qc_file)
        path(marker_file)
        path(tf_file)
        path(cyto_file)
    output:
        path("${id}*.pdf", emit: plots)
        path("${id}*.tsv", emit: summary)
    script:
        cmd = 
        """
        create_cnmf_summary.r \
        --spectra "${spectra_file}" \
        --output ${id}_${k} \
        --threshold ${params.cnmf.summarize.threshold}\
        """
        
        if (params.cnmf.summarize.topn != null) {
            cmd += " --topn ${params.cnmf.summarize.topn}"
        }
        
        if (params.cnmf.summarize.databases != null) {
            cmd += " --databases ${params.cnmf.summarize.databases}"
        }
        
        if (params.cnmf.summarize.tests != null) {
            cmd += " --tests ${params.cnmf.summarize.tests}"
        }
        
        if (params.cnmf.summarize.scale_spectra) {
            cmd += " --scaleSpectra"
        }
    
        if (enrich_file.getFileName().toString() != "NO_ENRICH") {
            cmd += " --enrichment ${enrich_file}"
        }
        
        if (marker_file.getFileName().toString() != "NO_MARKER") {
            cmd += " --annot ${marker_file}"
        }
        
        if (qc_file.getFileName().toString() != "NO_QC") {
            cmd += " --qc ${qc_file}"
        }
        
        if (tf_file.getFileName().toString() != "NO_TF") {
            cmd += " --tfFile ${tf_file}"
        }
        
        if (cyto_file.getFileName().toString() != "NO_CYTO") {
            cmd += " --cytoFile ${cyto_file}"
        }
    
        cmd
}



process cnmf_qc_summary {
    
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_selection", mode: 'symlink'
    
    input:
        val(id)
        path(files)
    output:
        tuple path("${id}*.tsv"), path("${id}*.pdf"), emit: qc_report
    script:
    
    cmd = 
    """
    create_cnmf_qc_summary.r \
    --inputs ${files.join(",")} \
    --out ${id}
    """
    cmd
}
