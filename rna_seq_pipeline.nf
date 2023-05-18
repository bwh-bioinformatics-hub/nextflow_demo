#!/usr/bin/env nextflow

// Get Input data - could be fastq or bam files
csv_file = file(params.input, checkIfExists: true)
ch_input = extract_data(csv_file)


// genome inputs 
ch_fasta = Channel.value(file(params.fasta))
ch_gtf   = Channel.value(file(params.gtf))

(bam_input, fastq_input) = ch_input.into(2)

if(params.input_type == 'bam'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_built.into(3)
}else if(params.input_type == 'fastq'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_input.into(3)
}

def summary = [:]
if(params.star)               summary ['STAR indices']     = params.star

process FASTQC_RAW {
    publishDir (
        path: "${params.outdir}/Reports/FastQC",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(base), file(fastq) from fastqc_reads

    output:
    file("*.{html,zip}") into fastqc_raw

    script:
    """
    fastqc -q $fastq --threads 12
    """
}


/*
================================================================================
                            BUILD INDICES
================================================================================
*/
process STAR_INDEX {
        
        publishDir (
        path: "${params.outdir}/STAR/index",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file("STARIndex") into star_built

    script:
    """
    mkdir -p STARIndex

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN 32 \\
        --sjdbOverhang ${params.sjdbOverhang} \\
        --sjdbGTFfile $gtf \\
        --genomeDir STARIndex/ \\
        --genomeFastaFiles $fasta
    """
}

ch_star = params.star ? Channel.value(file(params.star)) : star_built

/*
================================================================================
                                    TRIM
================================================================================
*/
process TRIM_GALORE {
    
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.fq.gz",
        saveAs: { params.save_qc_intermediates ? "quality_control/trimgalore/${it}" : null }

    when:
    params.trim_fastq

    input:
    tuple val(base), file(fastq) from trimming_reads

    output:
    tuple val(base), file('*.fq.gz') into trim_reads_ch, trim_reads_ch_fastqc
    //file(*) into trim_results

    script:
        """
    trim_galore \\
        --cores 20 \\
        --gzip \\
        --fastqc \\
        --paired \\
        ${fastq[0]} \\
        ${fastq[1]}
    """
}


aligner_reads = params.trim_fastq ? trim_reads_ch : raw_reads

process FASTQC_TRIMMED {

    publishDir (
        path: "${params.outdir}/Reports/FastQC_trimmed",
        mode: 'copy',
        overwrite: 'true',
    )


    input:
    tuple val(base), file(fastq) from trim_reads_ch_fastqc

    output:
    file ("*.{html,zip}") into fastqc_trimmed

    script:
    """
    fastqc -q $fastq --threads ${task.cpus}
    """
}

/*
================================================================================
                                STAR Alignment
================================================================================
*/
process STAR{
    tag "${base}"
        publishDir (
        path: "${params.outdir}/STAR/",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(base), file(reads) from aligner_reads
    file(star_idx) from ch_star

    output:
    tuple val(base), file("${base}/${base}.Aligned.sortedByCoord.out.bam") into star_bam_files 

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p ${base}

    STAR \\
        --twopassMode Basic \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --genomeDir ${star_idx} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix ${base}/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outBAMsortingBinsN 150 \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN 20 \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

/*
================================================================================
                            Auxiliary functions
================================================================================
*/

// Check integer
def isValidInteger(value){
    value instanceof Integer
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}


// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "error:  Invalid CSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "error: Cannot find supplied FASTQ or BAM input file. If using input method CSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Read input files from input CSV
def extract_data(csvFile){
    Channel
        .fromPath(csvFile)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_keys = ["Sample_ID", "Read1", "Read2", "Bam"]
        if(!row.keySet().containsAll(expected_keys)) exit 1, "error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2', 'Bam'."

        checkNumberOfItem(row, 4)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)
        def bam = row.Bam.matches('NA') ? 'NA' : return_file(row.Bam)

        if(samples == '' || read1 == '' || read2 == '' || bam == '') exit 1, "error: a field does not contain any information. Please check your CSV file"
        if(read1.matches('NA') && read2.matches('NA') && bam.matches('NA')) exit 1, "error: A row in your CSV file appears to have missing information."
        if( !read1.matches('NA') && !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
        if( !read2.matches('NA') && !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
        if( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "error: A specified BAM file either has a non-recognizable extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

        // output tuple mimicking fromFilePairs if FASTQ provided, else tuple for BAM
        if(bam.matches('NA')){
            if(read2.matches('NA')){
                [ samples, read1 ]
            }else{
                [ samples, [read1, read2] ]
            }
        }else{
            [ samples, bam ]
        }

        }
}