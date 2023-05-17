### QC

1. Quality statistics of raw sequence data

   The data volume and sequence quality of each sample sequence data (fastq) were counted with `fastqc` (v0.11.8).

   The statistics of raw sequence data of each sample are shown in the table below：

   {}

   Column name explanation:

   (1) Sample: Name of sequencing sample;
   (2) Length(bp): Average length of reads;
   (3) Reads: Number of sequence reads;
   (4)  GC%: Calculate the percentage of the total number of bases G and C in the total number of bases;
   (5) Q20(%): Calculate the percentage of bases with phred value greater than 20 in the total base;
   (6) Q30(%): Calculate the percentage of bases with phred value greater than 30 in the total base;

2. Quality control of sequence data

   In high-throughput sequencing, there are usually some sequencing errors such as point mutation, rRNA contamination, and the quality of the end of the sequence is relatively low. In order to obtain higher quality and more accurate biological information analysis results, it is necessary to optimize the original sequencing data.

   Analysis software：`bowtie2` (v2.4.4)    `fastp` (v0.23.1)

   Optimization steps and parameters:

   1. Remove rRNA contamination

      The reference index of `bowtie2` was constructed using the fasta sequence of rRNA of NCBI corresponding species, and then the rRNA sequence in analysis data was removed using `bowtie2` alignment.

   2. Remove low quality bases and reads

      `Fastp`,  a effective tool to integrate multiple sequence data quality control methods, can filter low-quality, too short and too many N reads,  and remove joint pollution and cut off sub sequences with low mean value and so on. In RSAP, we use the default parameters recommended by `fastp` for quality control.

   After quality control of sequence data of each sample, the statistics are as follows:

   {}

   Column name explanation:

   (1) Sample: Name of sequencing sample;
   (2) Reads: Number of sequence reads;
   
   (3) Bases: Total number of bases;

   (4) Length: Average length of reads;
   (5) Q20: Calculate the percentage of bases with phred value greater than 20 in the total base;
   (6) Q30: Calculate the percentage of bases with phred value greater than 30 in the total base;
   
   (7) GC: Calculate the percentage of the total number of bases G and C in the total number of bases;

   (8) ReadPassRate: Proportion of reads that pass quality control in the total number of original reads;