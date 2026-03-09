process genotypeGVCFs {

  label 'process_medium'
  tag "genotypeGVCFs"

  container 'broadinstitute/gatk:4.6.0.0'
  publishDir "${params.outdir}/VCF", mode: "copy"

  input:
    tuple path(cohort_gvcf), path(cohort_gvcf_tbi)
    path genome_files

  output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi"), emit: vcf

  script:
    """
    set -euo pipefail
    REF=\$(basename "${params.genome_file}")

    gatk GenotypeGVCFs \\
      -R "\$REF" \\
      -V ${cohort_gvcf} \\
      -O cohort.vcf.gz

    gatk IndexFeatureFile -I cohort.vcf.gz
    """
}