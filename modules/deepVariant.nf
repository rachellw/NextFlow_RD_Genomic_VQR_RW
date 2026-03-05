// Use newest nextflow dsl
nextflow.enable.dsl=2

process DEEPVARIANT {

   if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    
  tag "$sample_id"

  // DeepVariant official image (pin a version you like)
  container 'google/deepvariant:1.6.1'

  publishDir "${params.outdir}/DEEPVARIANT", mode: 'copy'

  input:
  tuple val(sample_id), path(bamFile), path(bamIndex)
  path indexFiles
  // optional: BED/intervals for WES, etc.
  path intervals optional true

  output:
  tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"),  emit: vcf
  tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: tbi
  // If you decide to run in gVCF mode you can emit gvcf as well (see note below)

  script:
  // DeepVariant model type: WGS, WES, PACBIO, ONT_R104, etc.
  def model = params.deepvariant_model ?: 'WGS'
  def threads = task.cpus ?: 12

  // In DV images, /opt/deepvariant/bin/run_deepvariant is the common entry
  // We write outputs into the task work dir (relative filenames!)
  def extraIntervalsArg = intervals ? "--regions ${intervals}" : ""
  def shards = task.cpus as int
  """
  set -euo pipefail

  echo "Running DeepVariant for ${sample_id}"
  echo "Model: ${model}"
  echo "Threads: ${threads}"

  /opt/deepvariant/bin/run_deepvariant \\
    --model_type=${model} \\
    --ref=${ref_fasta} \\
    --reads=${bam} \\
    --output_vcf=${sample_id}.deepvariant.vcf.gz \\
    --num_shards=${threads} \\
    ${extraIntervalsArg}

  echo "DeepVariant complete"
  """
}