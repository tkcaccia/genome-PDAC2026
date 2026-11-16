nextflow.enable.dsl=2

process PREPARE {
    cpus 1
    memory '2 GB'
    output:
    path 'prepare.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' prepare '${params.package_root}' '${params.run_root}' '${params.asset_root}'
    touch prepare.done
    """
    stub:
    'touch prepare.done'
}

process FASTQC {
    tag "$sample"
    cpus 4
    memory '8 GB'
    time '2h'
    input:
    tuple val(index), val(sample), path(ready)
    output:
    tuple val(index), val(sample), path('fastqc.done')
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 20_fastqc_2_5.sbatch ${index}
    touch fastqc.done
    """
    stub:
    'touch fastqc.done'
}

process ALIGN {
    tag "$sample"
    cpus 8
    memory '32 GB'
    time '8h'
    input:
    tuple val(index), val(sample), path(ready)
    output:
    tuple val(index), val(sample), path('alignment.done')
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 02_align_markdup_2_5.sbatch ${index}
    touch alignment.done
    """
    stub:
    'touch alignment.done'
}

process COVERAGE {
    tag "$sample"
    cpus 4
    memory '16 GB'
    time '2h'
    input:
    tuple val(index), val(sample), path(ready)
    output:
    path 'coverage.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 05_coverage_qc_2_5.sbatch ${index}
    touch coverage.done
    """
    stub:
    'touch coverage.done'
}

process CALL_VARIANTS {
    tag "$patient"
    cpus 8
    memory '32 GB'
    time '12h'
    input:
    tuple val(index), val(patient), val(alignment_barrier)
    output:
    tuple val(index), val(patient), path('calls.done')
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 03_deepsomatic_2_5.sbatch ${index}
    touch calls.done
    """
    stub:
    'touch calls.done'
}

process ANNOTATE {
    tag "$patient"
    cpus 4
    memory '16 GB'
    time '3h'
    input:
    tuple val(index), val(patient), path(ready)
    output:
    path 'annotation.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 04_vep_2_5.sbatch ${index}
    touch annotation.done
    """
    stub:
    'touch annotation.done'
}

process MSI_SCAN {
    cpus 1
    memory '12 GB'
    time '6h'
    input:
    path ready
    output:
    path 'scan.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 26_msi_scan.sbatch 1
    touch scan.done
    """
    stub:
    'touch scan.done'
}

process MSI_PAIRED {
    tag "$patient"
    cpus 4
    memory '16 GB'
    time '48h'
    maxForks 4
    input:
    tuple val(index), val(patient), val(alignment_barrier), path(scan_ready)
    output:
    path 'msi.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' batch '${params.run_root}' 27_msi_paired.sbatch ${index}
    touch msi.done
    """
    stub:
    'touch msi.done'
}

process MSI_SUMMARY {
    cpus 1
    memory '4 GB'
    input:
    val paired_barrier
    output:
    path 'msi_summary.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' msi_summary '${params.run_root}'
    touch msi_summary.done
    """
    stub:
    'touch msi_summary.done'
}

process ANALYSES {
    cpus 8
    memory '48 GB'
    time '48h'
    input:
    val annotation_barrier
    output:
    path 'analyses.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' analyses '${params.run_root}'
    touch analyses.done
    """
    stub:
    'touch analyses.done'
}

process FIGURE_INVENTORY {
    cpus 2
    memory '16 GB'
    time '4h'
    input:
    path analyses_ready
    path msi_ready
    output:
    path 'figures.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' figures '${params.run_root}'
    touch figures.done
    """
    stub:
    'touch figures.done'
}

process EXPORT_FIGURES {
    cpus 1
    memory '8 GB'
    input:
    path figures_ready
    output:
    path 'exports.done'
    script:
    """
    ${params.python} '${projectDir}/stage.py' export '${params.run_root}'
    touch exports.done
    """
    stub:
    'touch exports.done'
}

process VERIFY_RESULTS {
    cpus 1
    memory '4 GB'
    publishDir "${params.reports_dir}", mode: 'copy'
    input:
    path exports_ready
    val qc_barrier
    output:
    path 'analysis_audit.json'
    script:
    """
    ${params.python} '${projectDir}/stage.py' audit '${params.run_root}' analysis_audit.json
    """
    stub:
    """
    echo '{"test":"stub-only","scientific_validation":false}' > analysis_audit.json
    """
}

process AUDIT_EXISTING {
    cpus 1
    memory '4 GB'
    publishDir "${params.reports_dir}", mode: 'copy'
    output:
    path 'analysis_audit.json'
    script:
    """
    ${params.python} '${projectDir}/stage.py' audit '${params.asset_root}' analysis_audit.json
    """
}

workflow {
    if (params.mode == 'audit') {
        AUDIT_EXISTING()
    } else {
        if (params.mode != 'full') error "Use --mode full or --mode audit"
        if (!params.run_root) error "--run_root is required for a fresh full run"
        def manifest = file(params.manifest, checkIfExists: true)
        def rows = manifest.readLines().drop(1).findAll { it.trim() }.collect { it.split('\t') }
        def samples = rows.withIndex().collect { row, index -> [index + 1, row[0]] }
        def patients = rows.collect { it[1] }.unique().sort { it.toInteger() }.withIndex().collect { patient, index -> [index + 1, patient] }
        PREPARE()
        FASTQC(Channel.fromList(samples).combine(PREPARE.out))
        ALIGN(Channel.fromList(samples).combine(PREPARE.out))
        COVERAGE(ALIGN.out)
        def aligned = ALIGN.out.collect().map { 'all_alignments_complete' }
        CALL_VARIANTS(Channel.fromList(patients).combine(aligned))
        ANNOTATE(CALL_VARIANTS.out)
        MSI_SCAN(PREPARE.out)
        MSI_PAIRED(Channel.fromList(patients).combine(aligned).combine(MSI_SCAN.out))
        MSI_SUMMARY(MSI_PAIRED.out.collect())
        ANALYSES(ANNOTATE.out.collect())
        FIGURE_INVENTORY(ANALYSES.out, MSI_SUMMARY.out)
        EXPORT_FIGURES(FIGURE_INVENTORY.out)
        VERIFY_RESULTS(EXPORT_FIGURES.out, FASTQC.out.mix(COVERAGE.out).collect())
    }
}