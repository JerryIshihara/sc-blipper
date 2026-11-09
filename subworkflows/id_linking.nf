#!/usr/bin/env nextflow
include { fetch_gene_id_reference; invert_id_link } from "../processes/utils.nf"
include { filter_biotype } from "../processes/utils.nf"


workflow fetch_id_linker {

    take:
        version
        convert_params

    main:
    
        // Sanity check
        if (convert_params.output_namespace !in ['ensembl', 'gene_name']) {
            error "convert.output_namespace must be one of 'ensembl' or 'gene_name'"
        }    
        
        // Always fetch this even if not used

        if (convert_params.ensembl_file != null) {
            
            if (convert_params.id_linker == null) {
                error "If convert.ensembl_file is set, convert.id_linker must be set."
            }
            
            ensembl_file = file(convert_params.ensembl_file)
            if (!ensembl_file.exists()) {
                error "Supplied ensembl file does not exist: ${convert_params.ensembl_file}"
            }
            ensembl = Channel.value(ensembl_file)
            ensembl_reference = null
        } else {
            ensembl_reference = fetch_gene_id_reference(version)
            ensembl = ensembl_reference.ensembl
        }

        if (convert_params.id_linker) {
            id_file = file(convert_params.id_linker)
            if (!id_file.exists()) {
                error "Supplied id linker file does not exist: ${convert_params.id_linker}"
            }
            id_linker     = Channel.value(id_file)
            id_linker_inv = invert_id_link(id_linker)

        } else {
            if (convert_params.output_namespace == 'gene_name') {
                // ENSEMBL → gene name
                id_linker     = ensembl_reference.ensembl_to_name
                id_linker_inv = ensembl_reference.name_to_ensembl
            } else {
                // gene name → ENSEMBL
                id_linker     = ensembl_reference.name_to_ensembl
                id_linker_inv = ensembl_reference.ensembl_to_name
            }
        }
        
        if (convert_params.biotype_filter != null) {
            biotype_filtered = filter_biotype(ensembl)
            biotype_ensembl=biotype_filtered.biotype_ensembl
            biotype_gene_name=biotype_filtered.biotype_gene_name
        } else {
            biotype_ensembl=Channel.value(file("NO_SUBSET"))
            biotype_gene_name=Channel.value(file("NO_SUBSET"))
        }
        //}
        //else {
        //    // No mapping
        //    dummy_file = file("NO_MAPPING")
        //    id_linker     = Channel.value(dummy_file)
        //    id_linker_inv = Channel.value(dummy_file)
        // }

    emit:
        ensembl_reference = ensembl
        id_linker = id_linker
        id_linker_inv = id_linker_inv
        biotype_ensembl=biotype_ensembl
        biotype_gene_name=biotype_gene_name
}
