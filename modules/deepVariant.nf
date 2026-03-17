nextflow.enable.dsl=2

process DEEPVARIANT {

  label 'process_high'
  tag "$sample_id"

  container 'google/deepvariant:1.9.0'
  publishDir "${params.outdir}/DEEPVARIANT", mode: 'copy'

  input:
    tuple val(sample_id), path(bamFile), path(bamIndex)
    path indexFiles

output:
tuple val(sample_id),
      path("${sample_id}.g.vcf.gz"),
      path("${sample_id}.g.vcf.gz.tbi"),
      emit: gvcf

path("${sample_id}.vcf.gz"), optional: true, emit: vcf
path("${sample_id}.vcf.gz.tbi"), optional: true, emit: vcf_tbi

  script:
def model  = params.deepvariant_model ?: 'WGS'
def shards = (task.cpus ?: 8) as int

"""
set -euo pipefail

echo "Running DeepVariant (gVCF) for ${sample_id}"
echo "Model: ${model}"
echo "Shards: ${shards}"

REF=\$(basename "${params.genome_file}")
if [ ! -f "\$REF" ]; then
  REF=\$(find -L . -maxdepth 1 \\( -name "*.fasta" -o -name "*.fa" \\) | head -n 1)
fi

echo "Using reference: \$REF"
ls -lh

/opt/deepvariant/bin/run_deepvariant \\
  --model_type=${model} \\
  --ref="\$REF" \\
  --reads="${bamFile}" \\
  --output_vcf="${sample_id}.vcf.gz" \\
  --output_gvcf="${sample_id}.g.vcf.gz" \\
  --num_shards=${shards}

if [ ! -f "${sample_id}.g.vcf.gz.tbi" ]; then
  tabix -p vcf "${sample_id}.g.vcf.gz"
fi

if [ ! -f "${sample_id}.vcf.gz.tbi" ]; then
  tabix -p vcf "${sample_id}.vcf.gz"
fi

echo "DeepVariant complete"
"""
}
  