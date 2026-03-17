// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
          DNASeq Pipeline Configuration
    ============================================
    platform        : ${params.platform}
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index genome    : ${params.index_genome}
    qsr truth vcfs  : ${params.qsrVcfs}
    output directory: ${params.outdir}
    fastqc          : ${params.fastqc}
    aligner         : ${params.aligner}
    variant caller  : ${params.variant_caller}
    bqsr            : ${params.bqsr}
    fastp           : ${params.fastp}
    degraded_dna    : ${params.degraded_dna}
    variant_recalibration: ${params.variant_recalibration}
    identity_analysis: ${params.identity_analysis}
    fastq samplesheet        : ${params.samplesheet}
    bam no index samplesheet: ${params.bam_no_index_samplesheet}
    bam samplesheet         : ${params.bam_samplesheet}
    gvcf samplesheet        : ${params.gvcf_samplesheet}
    ============================================
""".stripIndent()

// Conditionally include modules
if (params.index_genome) {
    include { indexGenome } from './modules/indexGenome'
}

if (params.fastp) {
    include {FASTP } from './modules/readTrimming'
}

if (params.fastqc) {
    include { FASTQC } from './modules/FASTQC'
}
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
if (params.bqsr) {
    include { baseRecalibrator } from './modules/BQSR'
}
include { combineGVCFs } from './modules/processGVCFs'
include { genotypeGVCFs } from './modules/genotype_GVCF'
if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) {
    include { identityAnalysis } from './modules/identityAnalysis'
}
if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else {
    error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem' or 'bwa-aln'."
}
if (params.variant_caller == 'haplotype-caller') {
    include { haplotypeCaller } from './modules/haplotypeCaller'
 } else if  (params.variant_caller == 'deepvariant') {
    include { DEEPVARIANT } from './modules/deepVariant'  
} else {
    error "Unsupported variant caller: ${params.variant_caller}. Please specify 'haplotype-caller' or 'deepvariant'."
}

if (params.degraded_dna) {
    include { mapDamage2 } from './modules/mapDamage'
    include { indexMapDamageBam } from './modules/indexBam'
}



workflow FROM_FASTQ {

    log.info "Starting pipeline from FASTQ input"

    // User decides to index genome or not
    if (params.index_genome) {
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Create qsrc_vcf_ch channel
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }

    read_pairs_ch.view { x -> "READ_PAIRS_CH -> ${x}" }

    // Run FASTP on read pairs
    if (params.fastp) {
        FASTP(read_pairs_ch)
        read_pairs_ch = FASTP.out.trimmedreads
    }

    // Run FASTQC on read pairs
    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }

    // Align reads to the indexed genome
    if (params.aligner == 'bwa-mem') {
        align_ch = alignReadsBwaMem(read_pairs_ch, indexed_genome_ch.collect())
    }
    else if (params.aligner == 'bwa-aln') {
        align_ch = alignReadsBwaAln(read_pairs_ch, indexed_genome_ch.collect())
    }

    align_ch.view { x -> "ALIGN_CH -> ${x}" }

    // Sort BAM files
    sort_ch = sortBam(align_ch)

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files
    indexed_bam_ch = indexBam(mark_ch)

    indexed_bam_ch.view { x -> "INDEXED_BAM_CH -> ${x}" }

    // Conditionally run mapDamage if degraded_dna parameter is set
    if (params.degraded_dna) {
        pre_mapDamage_ch = mapDamage2(indexed_bam_ch, indexed_genome_ch.collect())
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    }
    else {
        mapDamage_ch = indexed_bam_ch
    }

    // Create known sites channel for BQSR
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { file ->
            file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf')
        }
        .map { file -> "--known-sites ${file.getName()}" }
        .collect()

    // Run or skip BQSR
    if (params.bqsr) {
        bqsr_ch = baseRecalibrator(
            mapDamage_ch,
            knownSites_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    }
    else {
        bqsr_ch = mapDamage_ch
    }

    // Variant caller
    if (params.variant_caller == 'haplotype-caller') {
        hc_result = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = hc_result.gvcf
    }
    else if (params.variant_caller == 'deepvariant') {
        dv_result = DEEPVARIANT(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = dv_result.gvcf
    }
    else {
        error "Unsupported variant caller: ${params.variant_caller}"
    }

    gvcf_ch.view { x -> "GVCF_CH -> ${x}" }

    // Combine sample gVCFs into cohort lists
    all_gvcf_ch = gvcf_ch
        .collect()
        .map { rows ->
            tuple(
                rows.collect { it[0] },
                rows.collect { it[1] },
                rows.collect { it[2] }
            )
        }

    all_gvcf_ch.view { x -> "ALL_GVCF_CH -> ${x}" }

    // Combine GVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, indexed_genome_ch.collect())
    combined_gvcf_ch.view()

    // Joint genotyping
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, indexed_genome_ch.collect())
    final_vcf_ch.view()

    // Conditionally apply variant recalibration or filtering
    if (params.variant_recalibration) {

        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels'      : 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38'                           : 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38'                        : 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38'    : 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138'          : 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38' : 'known=true,training=true,truth=true,prior=12.0'
        ]

        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file ->
                file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf')
            }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                return "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            final_vcf_ch,
            knownSitesArgs_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    }
    else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())
    }

    filtered_vcf_ch.view()

    // Optional identity analysis
    if (params.identity_analysis) {

        psam_info_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row ->
                if (row.size() == 4) {
                    tuple(row[0], row[3])
                } else if (row.size() == 3) {
                    tuple(row[0], row[2])
                } else {
                    error "Unexpected row format in samplesheet: $row"
                }
            }

        def combined_psam_content = new StringBuilder("#IID\tSID\tPAT\tMAT\tSEX\n")

        psam_file_ch = psam_info_ch.map { sample_info ->
            def sample_id = sample_info[0]
            def sex = sample_info[1]

            if (!sex) {
                sex = "NA"
            }

            def sample_line = "${sample_id}\t${sample_id}\t0\t0\t${sex}".stripIndent().trim()
            combined_psam_content.append(sample_line + "\n")
        }

        psam_file_ch.subscribe {
            def combined_psam_file = new File("/tmp/combined_samples.psam")
            combined_psam_file.text = combined_psam_content.toString()
            return combined_psam_file
        }

        identity_analysis_ch = identityAnalysis(filtered_vcf_ch, psam_file_ch)
    }
}

workflow FROM_DEDUP_BAM_NO_INDEX {

    log.info "Starting pipeline from deduplicated BAM (no index)"

    // Reference genome/index handling
    if (params.index_genome) {
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Resource VCFs
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Samplesheet format:
    // sample_id    bam
    bam_ch = Channel
        .fromPath(params.bam_no_index_samplesheet)
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            tuple(row.sample_id, file(row.bam))
        }

    bam_ch.view { x -> "DEDUP_BAM_CH -> ${x}" }

    // Index BAM first
    indexed_bam_ch = indexBam(bam_ch)

    indexed_bam_ch.view { x -> "INDEXED_BAM_CH -> ${x}" }

    // Optional degraded DNA branch
    if (params.degraded_dna) {
        pre_mapDamage_ch = mapDamage2(indexed_bam_ch, indexed_genome_ch.collect())
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    }
    else {
        mapDamage_ch = indexed_bam_ch
    }

    // Known sites for BQSR
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { f ->
            f.getName().endsWith('.vcf.gz') || f.getName().endsWith('.vcf')
        }
        .map { f -> "--known-sites ${f.getName()}" }
        .collect()

    // Optional BQSR
    if (params.bqsr) {
        bqsr_ch = baseRecalibrator(
            mapDamage_ch,
            knownSites_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    }
    else {
        bqsr_ch = mapDamage_ch
    }

    // Variant caller
    if (params.variant_caller == 'haplotype-caller') {
        hc_result = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = hc_result.gvcf
    }
    else if (params.variant_caller == 'deepvariant') {
        dv_result = DEEPVARIANT(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = dv_result.gvcf
    }
    else {
        error "Unsupported variant caller: ${params.variant_caller}"
    }

    gvcf_ch.view { x -> "GVCF_CH -> ${x}" }

    // Bundle sample gVCFs into cohort lists
    all_gvcf_ch = gvcf_ch
        .collect()
        .map { rows ->
            tuple(
                rows.collect { it[0] },
                rows.collect { it[1] },
                rows.collect { it[2] }
            )
        }

    all_gvcf_ch.view { x -> "ALL_GVCF_CH -> ${x}" }

    // Combine GVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, indexed_genome_ch.collect())
    combined_gvcf_ch.view()

    // Joint genotyping
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, indexed_genome_ch.collect())
    final_vcf_ch.view()

    // Final filtering / VQSR
    if (params.variant_recalibration) {

        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels'      : 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38'                           : 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38'                        : 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38'    : 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138'          : 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38' : 'known=true,training=true,truth=true,prior=12.0'
        ]

        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file ->
                file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf')
            }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                return "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            final_vcf_ch,
            knownSitesArgs_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    }
    else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())
    }

    filtered_vcf_ch.view()
}
workflow FROM_DEDUP_BAM {

    log.info "Starting pipeline from deduplicated BAM inputs"

    // Reference genome/index handling
    if (params.index_genome) {
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Resource VCFs
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Read BAM samplesheet: sample_id, bam, bai
    bam_ch = Channel
        .fromPath(params.bam_samplesheet)
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            if (row.size() < 3) {
                error "BAM samplesheet must contain: sample_id, bam, bai"
            }
            tuple(row[0], file(row[1]), file(row[2]))
        }

    bam_ch.view()

    // Optional degraded DNA branch
    if (params.degraded_dna) {
        pre_mapDamage_ch = mapDamage2(bam_ch, indexed_genome_ch.collect())
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    } else {
        mapDamage_ch = bam_ch
    }

    // Known sites for BQSR
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { f ->
            f.getName().endsWith('.vcf.gz') || f.getName().endsWith('.vcf')
        }
        .map { f -> "--known-sites ${f.getName()}" }
        .collect()

    // Optional BQSR
    if (params.bqsr) {
        bqsr_ch = baseRecalibrator(
            mapDamage_ch,
            knownSites_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    } else {
        bqsr_ch = mapDamage_ch
    }

    // Variant caller
    if (params.variant_caller == 'haplotype-caller') {
        hc_result = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = hc_result.gvcf
    }
    else if (params.variant_caller == 'deepvariant') {
        dv_result = DEEPVARIANT(bqsr_ch, indexed_genome_ch.collect())
        gvcf_ch = dv_result.gvcf
    }
    else {
        error "Unsupported variant caller: ${params.variant_caller}"
    }

    // Debug
    gvcf_ch.view { x -> "GVCF_CH -> ${x}" }
    gvcf_ch.view { x -> "GVCF_CH CLASS=${x.getClass().name} VALUE=${x}" }

    // Combine sample gVCFs into cohort lists
    all_gvcf_ch = gvcf_ch
    .collect()
    .map { rows ->
        rows.each { r ->
            assert r instanceof List : "gvcf_ch item is not a tuple/list: ${r} (${r.getClass().name})"
            assert r.size() == 3 : "gvcf_ch item does not have 3 elements: ${r}"
        }

        tuple(
            rows.collect { it[0] },
            rows.collect { it[1] },
            rows.collect { it[2] }
        )
    }

    // Combine GVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, indexed_genome_ch.collect())

    // Joint genotyping
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, indexed_genome_ch.collect())

    // Final filtering or VQSR
    if (params.variant_recalibration) {

        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels'       : 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38'                            : 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38'                         : 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38'     : 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138'           : 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38'  : 'known=true,training=true,truth=true,prior=12.0'
        ]

        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file ->
                file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf')
            }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                return "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            final_vcf_ch,
            knownSitesArgs_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )
    }
    else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())
    }

    filtered_vcf_ch.view()
}
workflow FROM_GVCF {

    log.info "Starting pipeline from gVCF inputs"

    // Reference genome/index
    if (params.index_genome) {
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    } else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Read gVCF samplesheets
   gvcf_ch = Channel
    .fromPath(params.gvcf_samplesheet)
    .splitText()
    .map { line -> line.trim() }
    .filter { line -> line }
    .map { line ->
        def row = line.split('\t')

        if (row.size() != 3) {
            error "gVCF samplesheet must have exactly 3 tab-separated columns: sample_id, gvcf, gvcf_index. Got: ${line}"
        }

        tuple(
            row[0].trim(),
            file(row[1].trim()),
            file(row[2].trim())
        )
    }

gvcf_ch.view { x -> "GVCF_CH -> class=${x.getClass().name} value=${x}" }
    // Combine sample gVCFs into cohort lists
    all_gvcf_ch = gvcf_ch
        .collect()
        .map { rows ->
            tuple(
                rows.collect { it[0] },
                rows.collect { it[1] },
                rows.collect { it[2] }
            )
        }

    all_gvcf_ch.view()

    // Combine gVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, indexed_genome_ch.collect())

    // Joint genotyping
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, indexed_genome_ch.collect())

    final_vcf_ch.view()

    // Variant filtering
    if (params.variant_recalibration) {

        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels'      : 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38'                           : 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38'                        : 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38'    : 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138'          : 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38' : 'known=true,training=true,truth=true,prior=12.0'
        ]

        qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file ->
                file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf')
            }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                return "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            final_vcf_ch,
            knownSitesArgs_ch,
            indexed_genome_ch.collect(),
            qsrc_vcf_ch.collect()
        )

    } else {

        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())

    }

    filtered_vcf_ch.view()
}
workflow FASTQC_only {
    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()

    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }
}

workflow {

    log.info "Selecting pipeline start point from provided inputs"

    def provided_inputs = [
        params.samplesheet,
        params.bam_no_index_samplesheet,
        params.bam_samplesheet,
        params.gvcf_samplesheet
    ].findAll { it != null && it.toString().trim() }

    if (provided_inputs.size() == 0) {
        error """
        No valid input supplied.
        Please provide one of:
          --samplesheet
          --bam_no_index_samplesheet
          --bam_samplesheet
          --gvcf_samplesheet
        """.stripIndent()
    }

    if (provided_inputs.size() > 1) {
        error """
        Multiple input types were supplied.
        Please provide only one of:
          --samplesheet
          --bam_no_index_samplesheet
          --bam_samplesheet
          --gvcf_samplesheet
        """.stripIndent()
    }

    if (params.gvcf_samplesheet) {
        log.info "Detected gVCF samplesheet -> starting FROM_GVCF"
        FROM_GVCF()
    }
    else if (params.bam_samplesheet) {
        log.info "Detected indexed deduplicated BAM samplesheet -> starting FROM_DEDUP_BAM"
        FROM_DEDUP_BAM()
    }
    else if (params.bam_no_index_samplesheet) {
        log.info "Detected deduplicated BAM samplesheet without index -> starting FROM_DEDUP_BAM_NO_INDEX"
        FROM_DEDUP_BAM_NO_INDEX()
    }
    else if (params.samplesheet) {
        log.info "Detected FASTQ samplesheet -> starting FROM_FASTQ"
        FROM_FASTQ()
    }
}
workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}
