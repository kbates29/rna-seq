#! /home/kameron/bin/nextflow

/*
 * Starting with a raw fastq performs adapter trimming and QC with Trim Galore, 
 * aligns to reference genome with STAR and obtains gene counts with featurecounts
 */


/*
 * Defines parameters
 */

params.reads_dir =  null
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
    .fromFilePairs(params.reads_dir/params.reads)
    .ifEmpty {error "Cannot find reads matching: ${params.reads}"}
    .into {reads_fast_qc; reads_trim}


/*
 * Run fastqc prior to trimming to compare afterwords
 */
process fast_qc{
    label 'low_memory'
    tag "$sampleId"
    publishDir "${params.out_dir}/fastqc", mode: 'copy', 
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$sampleId/$filename" : "$filename"}

    input: 
        set sampleId, file(reads) from reads_fast_qc
    
    output:
        file "*_fastqc.{zip,html}" into fastqc_results

    script:
        """
        fastqc -q $reads
        """
}



/*
 *  Run Trim Galore which trims adapters with cutadapt and then reruns fastqc
 */

process trim_galore{
    label 'low_memory'
    tag "$sampleId"
    publishDir "${params.out_dir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$sampleId/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$sampleId/$filename"
        }

    input:
        set sampleId, file(reads) from reads_trim

    output:
        file "*.gz" into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc
    
    script:
        """
        trim_galore --paired --fastqc -a $forward -a2 $reverse $reads
        """
}


/*
 * Run star 
 */
 process star{
    label 'star'
    tag "$sampleId"
    publishDir "${$params.out.dir}/star/$sampleId",

    input:
        set sampleId, file(reads) from trimmed_reads
        file gtf from params.gtf

    output:
        set file("*.Log.final.out"), file("*.bam") into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log
        file "${sampleId}.Aligned.sortedByCoord.out.bam.bai" into aligned_index

    script:
        """
        STAR --runThreadN 8 --genomeDir ${params.genome_dir} \
        --sjdbGTFfile ${params.genome_dir}/$gtf \
        --outFileNamePrefix $sampleId --outSamtype Bam SortedByCoordinate \
        --twopassMode Basic --readFilesCommand zcat
        """
 }


 process featurecounts{
    labal 'low_memory'
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

 process merge_counts{
    label 'low_memory'
    publishDir "${params.out.dur}/featurecounts"

    input:
        file input_counts from count_files.collect()
    
    output:
        file "merged_counts.txt"
 }