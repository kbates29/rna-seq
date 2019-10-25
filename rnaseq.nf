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
params.gff = null
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
 * Run star and index bam files
 */
 process star{
    executor "sge"
    clusterOptions="-S /bin/bash -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/star", mode: 'copy', 
        saveAs: { filename ->
            if (filename.indexOf("Log") >0 ) "logs/$filename"
            else if (filename.indexOf("_STARgenome") >0) "STARgenome/$filename"
            else if (filename.indexOf("_STARpass1") >0) "STARpass1/$filename"
            else ("$filename")
        }
       
    input:
        file reads from trimmed_reads

    output:
        file "*.bam" into star_aligned, create_config
        file "*SJ.out.tab" into splice_junction
        file "*.out" into alignment_logs
        file "*Log.out" into star_log
        file "*_STARgenome"
        file "*_STARpass1"
        file "*.bai" into index_bam
        

    script:
        sampleId = reads[0].toString() - ~/(_1)?(\.txt)?(\.gz)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.txt)?(\.gz)?$/
        sampleId = sampleId + "_"
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

        samtools index ${sampleId}Aligned.sortedByCoord.out.bam
        """
 }

/*
 *  Run featurecounts 
 */
 process featurecounts{
    executor "sge"
    clusterOptions="-S /bin/bash -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/counts", mode: 'copy'
 
    input:
        file bam from star_aligned

    output:
        file "*counts.txt" into sample_counts

    script:
        sampleId = bam[0].toString() - 'Aligned.sortedByCoord.out.bam'
        """
        export PATH=$PATH

        featureCounts -p -T 8 -t exon -a ${params.genome_dir}/${params.gtf}\
        -o ${sampleId}counts.txt $bam
        """
 }


/*
 * Merge counts into one file 
 */
 process mergecounts{
    executor "sge"
    clusterOptions="-S /bin/bash"
    publishDir "${params.out_dir}/counts", mode: 'copy'   
   
    input:
        file input_files from sample_counts.collect()
   
    output:
        file 'merged_counts.txt'

    script:
        gene_id = "<(tail -n +2 ${input_files[0]} | cut -f1)"
        counts = input_files.collect{filename ->
            "<(tail -n +2 ${filename} | cut -f7)"}.join("\t")
        """
        paste $gene_id \t $counts > merged_counts.txt
        """
 }

/*
 * create experiemnt info for majiq config
 */
 process experiment_config {
    executor "sge"
    clusterOptions="-S /bin/bash"

    input:
        file bam from create_config.collect()
   
    output:
        file "experiment.txt" into experiment
   
   script:

        """
        for file in ${bam}
        do
            key=\${file%%_*}
            value=\${file%.*}
            echo "\$key=\$value" >> experiment.txt
        done 
        """
 }




/*
 * Creates majiq config file
 */
 process majiq_config {
   executor "sge"
   clusterOptions="-S /bin/bash"
   publishDir "${params.out_dir}/majiq", mode: 'copy' 
    
    input:
       file info from experiment

    output:
       file "majiq.config" into majiq_config

    script:
       
       """
       printf '[info]\nreadlen=149\nbamdirs=${params.out_dir}star/\ngenome=GRCh38\nstrandness=None\n[experiments]\n' > majiq.config
       cat $info >>majiq.config
       """
 }

/*
 * Majiq build
 */
 process majiq_build {
    executor "sge"
    clusterOptions="-S /bin/bash -pe smp 32"
    publishDir "${params.out_dir}/majiq", mode: 'copy'
     
    input:
        file config from majiq_config
     
    output:
        file "*.majiq" into majiq mode flatten
        file "*.sj" into sj_out
        file "*splicegraph.sql" into splicegraph
     
    script:
    """
    export PATH=$PATH

    majiq build ${params.genome_dir}/${params.gff} -c $config -j 32 -o ./
    """
 }

/*
 * Majiq psi
 */
 process majiq_psi {
    executor "sge"
    clusterOptions="-S /bin/bash -pe smp 32"
    tag "$name"
    publishDir "${params.out_dir}/majiq/psi", mode: 'copy'

    input:
        file mj from majiq

    output:
        file "*.voila" into voila
        file "*.tsv" into tsv
        file "*.log" into logs
        
    
    script:
        name = mj[0].toString() -'_Aligned.sortedByCoord.out.majiq'
        """
        export PATH=$PATH

        majiq psi $mj -n $name -j 32 -o ./ 
        """
 }

