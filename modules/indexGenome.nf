/*
 * Define the indexGenome process that creates a BWA index
 * given the genome fasta file
 */
process indexGenome {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'variantvalidator/indexgenome:1.1.0'

    publishDir("${params.outdir}/GENOME_IDX", mode: "copy")

    input:
    path genomeFasta

    output:
    tuple path(genomeFasta),
          path("${genomeFasta}.*"),
          path("${genomeFasta.baseName}.dict")

    script:
    """
    set -euo pipefail

    echo "Running Index Genome"

    # Generate BWA index
    bwa index "${genomeFasta}"

    # Generate samtools faidx
    samtools faidx "${genomeFasta}"

    # Generate GATK/Picard sequence dictionary with correct name
    picard CreateSequenceDictionary \
        R="${genomeFasta}" \
        O="${genomeFasta.baseName}.dict"

    echo "Genome Indexing complete."
    """
}
