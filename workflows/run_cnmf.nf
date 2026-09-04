#!/usr/bin/env nextflow

// Processes
include { cnmf_prepare; cnmf_factorize; cnmf_combine; cnmf_kselection; cnmf_consensus; cnmf_ktree; cnmf_summarize; cnmf_qc_summary } from "../processes/cnmf.nf"
include { gsea; ora; decoupler; fetch_omnipath } from "../processes/enrichment.nf"
include { magma_assoc } from "../processes/magma.nf"
include { magma_concat as magma_concat_main } from "../processes/magma.nf"
include { magma_concat as magma_concat_per_k } from "../processes/magma.nf"
include { concat_enrichment_results as concat_enrichment_results_main } from "../processes/enrichment.nf"
include { concat_enrichment_results as concat_enrichment_results_per_k } from "../processes/enrichment.nf"

// Subworkflows
include { fetch_id_linker } from "../subworkflows/id_linking.nf"
include { magma_base } from "../subworkflows/magma.nf"
include { cnmf_stage } from "../subworkflows/cnmf_stage.nf"

// Main workflow
workflow cnmf {
    //TODO: make gene linker file ids unique
    main:
        //------------------------------------------------------------
        // Determine the number of cnmf workers
        int total_cnmf_runs = params.cnmf.n_iter * params.cnmf.k.split(",").size()
        int cnmf_workers = 0
        
        if (params.cnmf.n_workers == null || params.cnmf.n_workers < 0 || params.cnmf.n_workers == 1) {
            cnmf_workers = (total_cnmf_runs - 1).toInteger()
        } else if (params.cnmf.n_workers > 0 && params.cnmf.n_workers < 1) {
            // fractional: interpret as proportion of total runs
            cnmf_workers = (total_cnmf_runs * params.cnmf.n_workers.toDouble()).toInteger()
        } else {
            cnmf_workers = params.cnmf.n_workers.toInteger()
        }
        
        if (cnmf_workers >= total_cnmf_runs || cnmf_workers < 2) {
            throw new Exception("Number of workers (cnmf.n_workers) cannot exceed number of cnmf runs (cnmf.n_iter * cnmf.k.size()) or be smaller then 2")
        }

        //------------------------------------------------------------
        // Id linking and ensembl reference
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        
        // Add the output in the global space for convenience
        is_ensembl = params.convert.output_namespace == 'ensembl'
    
        //--------------------------------------------------------
        //                       cNMF
        //--------------------------------------------------------
        // Stage and preprocess (batch correct) the cNMF files
        cnmf_staged = cnmf_stage(fetch_id_linker.out.id_linker,
                   fetch_id_linker.out.biotype_ensembl,
                   fetch_id_linker.out.biotype_gene_name)
        
        // Prepare the cNMF input folder. The process selects CPU or GPU-aware
        // arguments and environment from params.cnmf_gpu.prepare.
        cnmf_prepared = cnmf_prepare(cnmf_staged.cnmf_in,
                                     cnmf_staged.cnmf_in_tpm,
                                     cnmf_staged.cnmf_in_hvg)
         
        // Create the channel with the n_workers for jobs to run in parallel
        cnmf_factorize_in = cnmf_prepared
         .map { v -> (0..(cnmf_workers-1)).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        // Perform the factorizations. The process selects its CPU or GPU
        // runtime and arguments from the resolved cNMF configuration.
        cnmf_factorize_out = cnmf_factorize(cnmf_factorize_in, cnmf_workers)
        
        // Wait for all the workers to complete, then:
        // Flatten all the output into a single path object
        cnmf_factorize_out = cnmf_factorize_out.collect()
        
        // Combine the output   
        cnmf_combine_out = cnmf_combine(cnmf_prepared, cnmf_factorize_out)     
 
        // Make the consensus factorizations
        cnmf_consensus_in = cnmf_prepared
         .map { v -> (params.cnmf.k.split(",")).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        // Create consensus with or without h5ad
        if (params.cnmf.save_h5ad) {
            cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out, cnmf_staged.cnmf_in.map{row -> row[1]})
        } else {
            cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out, file("NO_H5AD"))
        }
        // This is the end of the cnmf processing
        
        
        // Summarize QC metrics
        qc_summary_out = cnmf_qc_summary(params.rn_runname, cnmf_out.qc_metrics.map{row -> row[1]}.toSortedList())
    
        // TODO: Wrap the below in : if(params.cnmf.k.split(",").size() >= 2)
        // This is to enable running with a single k value
        // Make the k selection plot
        cnmf_kselection(cnmf_prepared, cnmf_combine_out)
        
        // Plot k-selection tree (how programs relate to eachother)
        if (params.cnmf.ktree_plot) {
            // Collect all (k, file) pairs into a list
            all_files = cnmf_out.spectra_score.map{row -> row[1]}.toSortedList()

            // Map to the "k=fileName" strings and join with space
            pairs_string = cnmf_out.spectra_score
            .flatMap { row -> "${row[0]}=${row[1].getName()}" }
            .toSortedList()
            .map { list -> list.join(' ') }

            // Create a single emission channel for the next process with all files and the string
            cnmf_ktree(all_files, pairs_string, params.rn_runname)
        }
        

        //--------------------------------------------------------
        // Optional post-proccessing of cnmf
        //--------------------------------------------------------

        //----------------------------------------------------------------------------------
        // Run gene set enrichment
        if (params.cnmf.run_enrichment) {
            
            // Channel with gmt files
            if (params.enrich.gmt_files != null) {
                if (params.enrich.gmt_files == 'DEFAULT') {
                    if (is_ensembl) {
                        gmt_files = Channel.fromPath("${projectDir}/assets/gene_sets/ensembl/*.gmt").collect()
                    } else {
                        gmt_files = Channel.fromPath("${projectDir}/assets/gene_sets/symbols/*.gmt").collect()
                    }
                } else {
                    gmt_files = Channel.from(params.enrich.gmt_files.split(",")).collect()
                }
            } else {
                if (params.cnmf.run_gsea || params.cnmf.run_ora) {
                    throw new Exception('No enrich.gmt files specified, but cnmf.run_gsea and or cnmf.run_ora is true')
                }
            }

            // Universe channel
            if (params.enrich.universe) {
                universe = file(params.enrich.universe)
                if (!universe.exists()) {
                    throw new Exception("Supplied id universe file does not exist")
                }
                // Custom conversion file
                universe = Channel.value(file(params.enrich.universe))
            } else {
                // In case there is no universe file
                universe = Channel.value(file("NO_UNIVERSE"))
            }
    
            if (params.cnmf.k_ignore != null) {
                gsea_in = cnmf_out.spectra_score.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                gsea_in = cnmf_out.spectra_score.map{i -> tuple("k_"+i[0], i[1])}
            }

            //----------------------------------------------------------------------------------
            // Run GSEA
            if (params.cnmf.run_gsea) {
                gsea_out = gsea("cnmf/consensus/${params.rn_runname}",
                    gsea_in,
                    gmt_files,
                    universe,
                    true)
            } else {
                gsea_out = [:]
                gsea_out.std = Channel.empty()
            }
            
            //----------------------------------------------------------------------------------
            // Run ORA
            if (params.cnmf.run_ora) {
                ora_out = ora("cnmf/consensus/${params.rn_runname}",
                    gsea_in,
                    gmt_files,
                    universe,
                    true,
                    false,
                    0,
                    false,
                    params.enrich.use_top)   
            } else {
                ora_out = [:]
                ora_out.std = Channel.empty()
            }

            //----------------------------------------------------------------------------------
            // Run decoupler (Progeny + collecTRI)
            if (params.cnmf.run_decoupler) {
                
                // Fetch the reference files        
                if (params.enrich.progeny_file != null && params.enrich.collectri_file != null) {
                    progeny_network = Channel.value(file(params.enrich.progeny_file))
                    collectri_network = Channel.value(file(params.enrich.collectri_file))
                } else if (params.enrich.progeny_file != null && params.enrich.collectri_file == null || params.enrich.progeny_file == null && params.enrich.collectri_file != null) {
                    error "Must specify both progeny and collectri files, or neither to fetch from omnipath"
                } else {
                    omnipath = fetch_omnipath()
                    progeny_network = omnipath.progeny
                    collectri_network = omnipath.collectri
                }
                
                // Filter which K values to ignore
                if (params.cnmf.k_ignore != null) {
                    decoupler_in = cnmf_out.spectra_score.filter { tuple ->
                        def (k, path) = tuple
                        !(k in params.cnmf.k_ignore.split(','))
                    }.map{i -> tuple("k_"+i[0], i[1])}
                } else {
                    decoupler_in = cnmf_out.spectra_score.map{i -> tuple("k_"+i[0], i[1])}
                }
                
                // Run decoupler with or without id conversion, depending on the input type
                if (is_ensembl) {
                    // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                    decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, fetch_id_linker.out.id_linker_inv, progeny_network, collectri_network) 
                } else {
                    // This assumes you have converted to gene symbols
                    decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, file("NO_MAPPING"), progeny_network, collectri_network)
                }
            } else {
                decoupler_out = [:]
                decoupler_out.std = Channel.empty()
            }

            
            //----------------------------------------------------------------------------------
            // Run magma
            if (params.cnmf.run_magma) {
                
                magma_base(params.magma, params.convert, params.rn_ensembl_version, fetch_id_linker.out.ensembl_reference)
                
                // Magma association with cnmf
                if (params.cnmf.k_ignore != null) {
                    magma_cnmf_in = cnmf_out.spectra_score.filter { tuple ->
                            def (k, path) = tuple
                            !(k in params.cnmf.k_ignore.split(','))
                        }.map{i -> tuple("k_"+i[0], i[1])}
                } else {
                    magma_cnmf_in = cnmf_out.spectra_score.map{i -> tuple("k_"+i[0], i[1])}
                }
                
                magma_assoc_in = magma_base.out.raw.combine(magma_cnmf_in)
                
                // Run regression based magma
                magma_assoc_out = magma_assoc(magma_assoc_in, true, universe, file("NO_MAPPING"))
                    
                // Concat the results in a single table
                concat_in = magma_assoc_out.out.collect().map{ list -> ["${params.rn_runname}", list]}
                magma_out = magma_concat_main("cnmf/consensus/", concat_in)
                
                // Collect the results per k value as well
                magma_out_per_k = magma_concat_per_k("cnmf/consensus/${params.rn_runname}", magma_assoc_out.per_database.groupTuple())

            } else {
                magma_out = [:]
                magma_out.std = Channel.empty()
                magma_out_per_k = [:]
                magma_out_per_k.std = Channel.empty()
            }

            //----------------------------------------------------------------------------------
            // Collect enrichment results for a summary
            // Merge all the output files and calculate FDR, optionally annotate

            // For cnmf workflow this doesn't make much sense to have, so default to not annotating
            // If annotation is needed, run cnmf, then manually run_enrich later
            annot_file = Channel.value(file("NO_ANNOT"))
            
            merge_in = Channel.empty()
            .mix(gsea_out.std.map{row -> row[1]},
                 ora_out.std.map{row -> row[1]},
                 decoupler_out.std.flatMap{row -> row[1]},
                 magma_out.std.map{row -> row[1]})
            .collect()
            .map{ list -> ["${params.rn_runname}", list]}
            
            concat_enrichment_results_main("cnmf/consensus/", merge_in, annot_file)

            // Per k value
            merge_in_per_k = Channel.empty()
            .mix(gsea_out.std,
                ora_out.std,
                decoupler_out.std.flatMap { id, files -> files.collect { file -> tuple(id, file) }},
                magma_out_per_k.std)
            .groupTuple(by:0)
            
            concat_enrichment_results_per_k("cnmf/consensus/${params.rn_runname}", merge_in_per_k, annot_file)
        }
        
        //----------------------------------------------------------------------------------
        // Make a sumamry table per value of k
        if (params.cnmf.run_enrichment) {
            summarize_in = cnmf_out.spectra_score.map{row -> tuple("k_${row[0]}", row[1])}.combine(concat_enrichment_results_per_k.out.std_nominal, by: 0)
        } else {
            summarize_in = cnmf_out.spectra_score.map{row -> tuple("k_${row[0]}", row[1], file("NO_ENRICH"))}
        }
        
        // Add the QC file
        summarize_in = summarize_in.combine(cnmf_out.qc_metrics.map{row -> tuple("k_${row[0]}", row[1])}, by: 0)

        // Compact selection of marker / tf / cyto files
        def chooseFile = { paramVal, defaultPath, noTag ->
            if (paramVal == null) Channel.value(file("NO_${noTag}"))
            else if (paramVal == "DEFAULT") Channel.value(file(defaultPath))
            else Channel.value(file(paramVal))
        }

        def suffix = is_ensembl ? "_ensembl.tsv" : "_gene_name.tsv"

        marker_file = chooseFile(params.cnmf.summarize.marker_file,
                                "${projectDir}/assets/markers/CD4_markers${suffix}",
                                "MARKER")

        tf_file = chooseFile(params.cnmf.summarize.tf_file,
                            "${projectDir}/assets/markers/lambert_2018_tfs${suffix}",
                            "TF")

        cyto_file = chooseFile(params.cnmf.summarize.cyto_file,
                            "${projectDir}/assets/markers/cytokines${suffix}",
                            "CYTO")
                       
        cnmf_summarize(params.rn_runname, summarize_in, marker_file, tf_file, cyto_file)
}
