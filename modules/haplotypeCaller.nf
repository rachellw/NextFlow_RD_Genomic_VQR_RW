nextflow.enable.dsl = 2

process haplotypeCaller {

    label 'process_low'
    container 'variantvalidator/gatk4:4.3.0.0'
    tag "$sample_id"
    publishDir "${params.outdir}/gvcf", mode: 'copy'
    input:
    tuple val(sample_id), path(bamFile), path(bamIndex)
    path indexFiles

    output:
    tuple val(sample_id),
          path("${sample_id}.g.vcf.gz"),
          path("${sample_id}.g.vcf.gz.tbi"),
          emit: gvcf

    script:
    """
    set -euo pipefail

    echo "Running HaplotypeCaller for sample: ${sample_id}"

    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename "${params.genome_file}")
    else
        genomeFasta=\$(find -L . -maxdepth 1 -name '*.fasta' -o -name '*.fa' | head -n 1)
    fi

    echo "Genome file: \${genomeFasta}"

    # Rename the dictionary file to the expected name if needed
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict" || true
    fi

    outputVcf="${sample_id}.g.vcf.gz"

    gatk HaplotypeCaller \\
        -R "\${genomeFasta}" \\
        -I ${bamFile} \\
        -O "\${outputVcf}" \\
        -ERC GVCF \\
        -A BaseQuality \\
        -A DepthPerSampleHC \\
        -A MappingQuality \\
        -A QualByDepth \\
        -A MappingQualityRankSumTest \\
        -A ReadPosRankSumTest \\
        -A FisherStrand \\
        -A StrandOddsRatio \\
        -A MappingQualityZero \\
        -A InbreedingCoeff \\
        -A BaseQualityRankSumTest \\
        -A HaplotypeFilteringAnnotation

    echo "Sample: ${sample_id} VCF: \${outputVcf}"

    gatk IndexFeatureFile -I "\${outputVcf}"
    
    echo "Variant Calling for sample: ${sample_id} complete"
    """
}
