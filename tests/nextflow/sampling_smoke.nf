// Exercise production processes with synthetic inputs; no farm or annotation downloads.
nextflow.enable.dsl = 2
include { cnmf_prepare; cnmf_sample; cnmf_factorize; cnmf_combine; cnmf_consensus; cnmf_kselection } from '../../processes/cnmf.nf'

workflow {
    input = Channel.value(tuple('smoke', file(params.test_counts), file(params.test_tpm), file(params.test_hvg)))
    if (params.cnmf.sampling.enabled) {
        samples = cnmf_sample(input)
        input = samples.map { id, folder ->
            def manifest = new groovy.json.JsonSlurper().parse(folder.resolve('manifest.json').toFile())
            def draw = manifest.samples[0]
            tuple(id, folder.resolve(draw.counts), folder.resolve(draw.tpm), folder.resolve(manifest.shared_hvg))
        }.first()
    }
    counts = input.map { id, counts, tpm, hvg -> tuple(id, counts) }
    fits = cnmf_prepare(counts, input.map { row -> row[2] }, input.map { row -> row[3] })
    workers = fits.flatMap { id, folder -> (0..1).collect { worker -> tuple(id, folder, worker) } }
    factors = cnmf_factorize(workers, 2).collect()
    merged = cnmf_combine(fits, factors)
    cnmf_kselection(fits, merged)
    consensus = fits.flatMap { id, folder -> params.cnmf.k.split(',').collect { k -> tuple(id, folder, k) } }
    h5ad = params.cnmf.save_h5ad ? input.map { row -> row[1] } : file(params.test_no_h5ad)
    cnmf_consensus(consensus, merged, h5ad)
}
