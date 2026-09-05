#!/usr/bin/env nextflow
/*
 * imputation-calibration-benchmark :: Phase 3 imputation workflow
 * Dual-condition Beagle 5.5 imputation of array-QC'd targets against
 * an ancestry-matched (A) and a deliberately mismatched (B) reference panel.
 * Module conventions follow biomarker-concordance-pipeline (meta map, versions).
 */
nextflow.enable.dsl = 2

params.gt         = "data/derived/array_targets.qc.vcf.gz"
params.panel_a    = "data/derived/panelA.vcf.gz"
params.panel_b    = "data/derived/panelB.vcf.gz"
params.map        = "data/raw/plink.chr20.GRCh37.map"
params.beagle_jar = "data/raw/beagle.27Feb25.75f.jar"
params.outdir     = "results"
params.xmx_gb     = 12
params.nthreads   = 6

process BEAGLE_IMPUTE {
    tag "${meta.id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(gt), path(ref), path(map), path(jar)

    output:
    tuple val(meta), path("imputed_${meta.id}.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    java -Xmx${params.xmx_gb}g -jar ${jar} \\
        gt=${gt} ref=${ref} map=${map} \\
        out=imputed_${meta.id} nthreads=${params.nthreads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(java -jar ${jar} 2>&1 | head -1 | sed 's/.*(version //;s/).*//')
        java: \$(java -version 2>&1 | head -1)
    END_VERSIONS
    """
}

workflow {
    conditions = Channel.of(
        [[id: 'A', condition: 'matched'],    file(params.gt), file(params.panel_a)],
        [[id: 'B', condition: 'mismatched'], file(params.gt), file(params.panel_b)],
    ).map { meta, gt, ref -> [meta, gt, ref, file(params.map), file(params.beagle_jar)] }

    BEAGLE_IMPUTE(conditions)
}
