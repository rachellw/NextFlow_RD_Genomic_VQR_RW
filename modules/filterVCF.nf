process filterVCF {

  tag "filterVCF"
  publishDir "${params.outdir}/FILTERED_VCF", mode: "copy"

  input:
  tuple path(vcf), path(vcf_index)
  path genomeFiles

  output:
  tuple path("filtered.vcf.gz"), path("filtered.vcf.gz.tbi")

  script:
  """
  REF=\$(basename ${params.genome_file})

  gatk VariantFiltration \
      -R \$REF \
      -V ${vcf} \
      -O filtered.vcf.gz

  gatk IndexFeatureFile -I filtered.vcf.gz
  """
}