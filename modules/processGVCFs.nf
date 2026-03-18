process combineGVCFs {

  label 'process_medium'
  tag "combineGVCFs"

  container 'broadinstitute/gatk:4.6.0.0'
  publishDir "${params.outdir}/GVCF", mode: "copy"

  input:
    tuple val(sample_ids), path(gvcfs), path(gvcf_indices)
    path genome_files

  output:
    tuple path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi"), emit: combined

  script:
    // Build -V args for each gVCF
    def v_args = gvcfs.collect { "-V ${it}" }.join(' ')

    """
    set -euo pipefail
    REF=\$(basename "${params.genome_file}")
    # GATK expects basename.dict, not basename.fasta.dict
    if [[ -f "${REF}.dict" && ! -f "${REF%.*}.dict" ]]; then
    ln -s "${REF}.dict" "${REF%.*}.dict"
    fi
    
    gatk CombineGVCFs \\
      -R "\$REF" \\
      ${v_args} \\
      -O cohort.g.vcf.gz

    # Ensure tabix index exists
    gatk IndexFeatureFile -I cohort.g.vcf.gz
    """
}