#!/usr/bin/env nextflow
include { convert_and_merge } from "../subworkflows/convert_merge.nf"
include { cnmf_sampling; cnmf_pre_process } from "../processes/cnmf.nf"
include { scvi_cnmf } from "../processes/scvi.nf"

workflow cnmf_stage {
    
    take:
        id_linker
        biotype_ensembl
        biotype_gene_name
    main:

        //------------------------------------------------------------ 
        // Prepare cnmf channels
        if (params.cnmf.input.h5ad != null && params.cnmf.input.tpm != null && params.cnmf.input.hvg != null) {
            // Override with previous output
            cnmf_in = Channel.value(tuple(params.rn_runname, file(params.cnmf.input.h5ad)))
            cnmf_in_tpm = Channel.value(file(params.cnmf.input.tpm))
            cnmf_in_hvg = Channel.value(file(params.cnmf.input.hvg))
        } else {
            //------------------------------------------------------------ 
            // Merge seurat and H5 files
            convert_and_merge(params.rn_manifest,
                                params.rn_runname,
                                id_linker,
                                params.convert.subset_genes,
                                biotype_ensembl,
                                biotype_gene_name,
                                params.convert.output_namespace)
                                
            merge_out = convert_and_merge.out.merge_out

            // Sample raw counts once, before preprocessing/scVI and HVG selection.
            if (params.sampling.enabled == true) {
                cnmf_sampling(merge_out)
                merge_out = cnmf_sampling.out.sampled.first()
            }
            
            //--------------------------------------------------------
            // Optionally run cNMF preprocessing, which creates harmony/scvi corrected data
            if (params.cnmf.preprocess) {
                // Preprocess the cnmf file
                if (params.preprocess.batch_correction == "scvi") {
                    cnmf_preprocess = scvi_cnmf(merge_out)
                } else {
                    cnmf_preprocess = cnmf_pre_process(merge_out)
                }
                
                // Channel with run name and .h5ad file with counts
                cnmf_in = cnmf_preprocess.preprocessed
                cnmf_in_tpm = cnmf_preprocess.normalized
                cnmf_in_hvg = cnmf_preprocess.hvg
            } else {    
                // Otherwise use the merged file, in case there is one batch,
                // just uses the other file
                cnmf_in = merge_out
                cnmf_in_tpm = Channel.value(file("NO_TPM"))
                cnmf_in_hvg = Channel.value(file("NO_HVG"))
            }
        }
        
    emit:
        cnmf_in=cnmf_in
        cnmf_in_tpm=cnmf_in_tpm
        cnmf_in_hvg=cnmf_in_hvg
}