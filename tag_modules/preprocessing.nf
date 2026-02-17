process PREPROCESS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq_dir)

    output:
    tuple val(sample_id), path("*.clean.fastq.gz")

    script:
    """
    tag_preprocessing.py -p -sample ${sample} -dir $fastq_dir -outdir $output_dir -primer $SA2 -barcode $barcode
    """
}

process ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    # Example: bwa (C-based tool)
    bwa mem -t ${task.cpus} ref.fa ${reads} | \
      samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    """
}

include { PREPROCESS } from './modules/preprocess'
include { ALIGN }      from './modules/align'
include { BAM_PROCESS } from './modules/bam_process'

workflow {

    Channel
        .fromFilePairs("data/*_{R1,R2}.fastq.gz", flat: true)
        .map { sample_id, reads -> tuple(sample_id, reads) }
        | PREPROCESS
        | ALIGN
        | BAM_PROCESS
}

