# rna-seq
Nexflow pipeline for RNA-seq experiment, currently implemented for use with SGE scheduler on HPC.
Aligns reads to reference genome, obtains gene counts and perorms alternative splicing analysis.

Step 1: Performs QC with FASTQC and Trim_galore

Step 2: Aligns reads using Star and sorts with Samtools

Step 3: Obtains read counts for each sample using Featurecounts

Step 4: Merges reads counts for all samples for downstream analysis

Step 5 & 6: Creates config file for Majiq

Step 7: Performs Majiq build.

Step 8: Performs splicing analysis per sample using Majiq
