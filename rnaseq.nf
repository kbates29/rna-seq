#!/gpfs_fs/home/batesk2/bin/nextflow

/*
 * Starting with a raw fastq performs adapter trimming and QC with Trim Galore, 
 * aligns to reference genome with STAR and obtains gene counts with featurecounts
 */


/*
 * Defines parameters
 */

params.reads = null
params.genome_dir = null
params.gtf = null
params.out_dir = null
params.forward_adap = null
params.reverse_adap = null

forward = params.forward_adap
reverse = params.reverse_adap

/*
 * Create channel for reads
 */
 Channel
    .fromFilePairs(params.reads)
    .ifEmpty {error "Cannot find reads matching: ${params.reads}"}
    .into {reads_fast_qc; reads_trim}


/*
 * Run fastqc prior to trimming to compare afterwords
 */
process fast_qc{
    executor 'sge'
    clusterOptions '-S /bin/bash'
    beforeScript 'printf "[%s] USER:\${USER:-none} JOB_ID:\${JOB_ID:-none} JOB_NAME:\${JOB_NAME:-none} HOSTNAME:\${HOSTNAME:-none} PWD:\$PWD\n"'
    tag "$sampleId"
    publishDir "${params.out_dir}/fastqc", mode: 'copy', 
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$sampleId/$filename" : "$filename"}

    input: 
        set sampleId, file(reads) from reads_fast_qc
    
    output:
        file "*_fastqc.{zip,html}" into fastqc_results

    script:
        """
        export PATH=$PATH
	    fastqc -q $reads
        """
}

/*
 *  Run Trim Galore which trims adapters with cutadapt and then reruns fastqc
 */

process trim_galore{
    executor "sge"
    clusterOptions="-S /bin/bash"
    beforeScript 'printf "[%s] USER:\${USER:-none} JOB_ID:\${JOB_ID:-none} JOB_NAME:\${JOB_NAME:-none} HOSTNAME:\${HOSTNAME:-none} PWD:\$PWD\n"'
    tag "$sampleId"
    publishDir "${params.out_dir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$sampleId/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$sampleId/$filename"
        }

    input:
        set val(sampleId), file(reads) from reads_trim

    output:
        file "*.gz" into trimmed_reads, test_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc
    
    script:
        """
        export PATH=$PATH
	    trim_galore --paired --fastqc -a $forward -a2 $reverse $reads
        """
}


/*
 * Run star 
 */
 process star{
    executor "sge"
    clusterOptions="-S /bin/bash -pe smp 8"
    beforeScript 'printf "[%s] USER:\${USER:-none} JOB_ID:\${JOB_ID:-none} JOB_NAME:\${JOB_NAME:-none} HOSTNAME:\${HOSTNAME:-none} PWD:\$PWD\n"'
    tag "$sampleId"
    publishDir "${params.out_dir}/star", mode: 'copy',
        saveAs: {filename ->
        if (filename.indexOf("Log.out")>0) "logs/$sampleId/$filename"
        }

    input:
        file reads from trimmed_reads

    output:
        file "*.bam" into star_aligned
        file "*Log.out" into star_log
        file "*.out" into alignment_logs
        file "*SJ.out.tab" into splice_junction
        

    script:
        sampleId = reads[0].toString() - ~/(_1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        export PATH=$PATH
        STAR --runMode alignReads --runThreadN 8 \
        --genomeDir ${params.genome_dir} \
        --readFilesIn $reads \
        --sjdbGTFfile ${params.genome_dir}/${params.gtf} \
        --sjdbOverhang 149 \
        --outFileNamePrefix $sampleId \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic --readFilesCommand zcat
        """
 }


 process featurecounts{
    labal 'low_mem'
    tag "$sampleId"
    publishDir "${params.out.dir}/featurecounts"

    input:
        file bam from star_aligned
        file gtf from params.gtf
    
    output:
        file "${sampleId}_featurecounts.txt" into count_files

    script:
        """
        featureCounts -T  -p -t exon -g gene_id -a $gtf -o ${sampleId}_featurecounts.txt $bam
        """
 }

/* process merge_counts{
    label 'low_memory'
    publishDir "${params.out.dur}/featurecounts"

    input:
        file input_counts from count_files.collect()
    
    output:
        file "merged_counts.txt"
 }*/
