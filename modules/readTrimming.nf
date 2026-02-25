/*
 * Run fastp on the [.  ] files
 */
nextflow.enable.dsl=2
process FASTP {

    label 'process_single'

    container 'staphb/fastp:1.1.0'

    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for the FASTP results
    publishDir("$params.outdir/FASTP", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*.fq.gz") , emit : trimmedreads
    path "*.json" , emit : json
    path "*.html" , emit : html 

    
    script:
    """
    echo "Running FASTP"
    mkdir -p fastp_${sample_id}_TR
    # Check the number of files in reads and run fastp accordingly
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
    fastp -i ${reads[0]} -I ${reads[1]} -o fastp_${sample_id}.R1.fq.gz -O fastp_${sample_id}.R2.fq.gz
    elif [ -f "${reads[0]}" ]; then
       -i fastp ${reads[0]} -o fastp_${sample_id}
    else
        echo "No valid [read files] found for sample ${sample_id}"
        exit 1
    fi

    echo "FASTP Complete"
    """
}
