# **TE Alignment Pipeline**

Pipeline to align RNAseq data and obtain counts for both TEs and genes

## **Outline**

* [Requirements](#requirements)
* [Pipeline](#pipeline)
  * [STEP 1 - Align uniquely mapped reads](#step-1---align-uniquely-mapped-reads)
  * [STEP 2 - Align mapped mapped reads](#step-2---align-multi-mapped-reads)
  * [STEP 3 (Optional) - Align uniquely mapped reads](#step-3-optional---create-a-log-for-all-the-alignments)
  * [STEP 4 - Create gtf files](#step-4---create-gtf-files)
    * [STEP 4.1 - Get the counts at the individual level](#step-41---get-the-counts-of-the-genes-and-tes-at-the-individual-level)
    * [STEP 4.2 - Get the counts at the subfamily level](#step-42---get-the-counts-of-the-genes-and-tes-at-the-subfamily-level)
* [Just Press Enter](#just-press-enter)
* [STAR Arguments explained](#star-arguments-explained)
* [featureCounts Arguments explained](#featurecounts-arguments-explained)

## **Requirements**

* STAR
* featureCounts
* a .gtf file with TEs annotated individually (you can see [here](#step-4---create-gtf-files) how to get it)
* a .gtf file with TEs annotated per subfamily (you can see [here](#step-4---create-gtf-files) how to get it)

## **Pipeline**

### STEP 1 - Align uniquely mapped reads

**Goal:** Get the genes and individually mapped TEs

```bash
# run STAR on paired-end FASTQ
STAR --runThreadN 16  \
     --genomeDir  [STAR_genome_index] \ # genome index
	 --limitSjdbInsertNsj 6000000 \
	 --alignIntronMax 500000 \
	 --alignMatesGapMax 500000 \
	 --alignEndsType EndToEnd \
     --outFileNamePrefix ${output_file}/${filename%_*}_ \
     --outFilterMultimapNmax 1 \
	 --outFilterMismatchNmax 999 \
	 --outFilterMismatchNoverLmax 0.03 \
	 --seedMultimapNmax 20000 \
	 --outSAMattributes All \
	 --outSAMmultNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn ${file%_*}_1.fastq ${file%_*}_2.fastq # For paired-end
```

Alternative, you could also `star_align_unique_TE.sh` (available in this repo)

```bash
bash star_align_unique_TE_passive.sh [folder_fastqs] [output_folder]
```

The file has code to run single or paired-end alignments, but it's not automatized, so you have to comment and uncomment which part you want. To run the in passive also uncomment the lines 2 and 3 :)

### STEP 2 - Align multi mapped reads

**Goal:** Get the subfamilies of TEs, to have a general idea of the rest of the TEs

It's almost the same arguments but with different values

```bash
# run STAR on paired-end FASTQ
STAR --runThreadN 72 \
	 --genomeDir  [STAR_genome_index] \ # genome index
	 --limitSjdbInsertNsj 6000000 \
	 --outFileNamePrefix ${output_file}/${filename%_*}_ \
	 --outFilterMultimapNmax 5000 \
	 --outFilterMismatchNmax 999 \
	 --outFilterMismatchNoverLmax 0.06 \
	 --seedMultimapNmax 20000 \
	 --outSAMattributes All \
	 --outSAMmultNmax 1 \
	 --outSAMtype BAM SortedByCoordinate \
	 --readFilesIn ${file%_*}_1.fastq ${file%_*}_2.fastq # For paired-end
```

Alternative, you could also run `star_align_multi_TE.sh` (available in this repo)

```bash
bash star_align_multi_TE.sh [folder_fastqs] [output_folder]
```

The file has code to run single or paired-end alignments, but it's not automatized, so you have to comment and uncomment which part you want. To run the in passive also uncomment the lines 2 and 3 :)

### STEP 3 (Optional) - Create a Log for all the alignments

If you want to track how many reads were unique/multi/unmapped for a batch of samples you can run the script `create_star_qc_log.sh`. This will give you a .tsv with the number of reads per sample in the STAR output folder

```bash
bash create_star_qc_log.sh [previous_STAR_output_folder]
```

### STEP 4 - Create gtf files

To create the custom made gtf I merged the normal annotation gtf (Gencode) with a repeatmasker with TEs

To get TEs file, I ran the `makeGTF.R` script from Miguel Casanova (https://github.com/milcs40/TEScripts), using as input the `hg38.fa.out` from repeat masker (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/)

According to him this is what the sript does:

> "The efficient quantification of the number of reads assigned to TEs depends on the accurate annotation of these elements across the genome and the creation of GTF files that allow counting the number of reads mapped to each annotated TE. To this end, we have developed a R function, makeGTF, which takes TE annotations from either UCSC or Repeatmasker databases, removes low complexity repeats (e.g. rRNAs, tRNAs, snRNAs, etc.), labels each individual TE annotation with an unique identifier and creates an attribute column that includes, besides the unique identifier, information about the subfamily, family and class to which the TE belongs. In addition, it adds the prefix “TE_” to each TE, allowing to easily distinguish between genes and TEs in downstream analysis. By default, the function generates three different types of GTF files: an unmodified GTF file, a GTF file that allows counting TE subfamilies (Table 4.1A) and a GTF file that allows counting individual TE instances (Table 4.1B). The GTF files are then used by featureCounts to scan the BAM files containing the genomic coordinates of each read aligned to the reference genome, and counts the number of reads mapping to TEs. The GTF files allow using either the subfamily or the unique TE identifier attributes to count the number of reads mapping to each subfamily or to individual TEs, respectively."

Running this will output:

* `rmsk_TEClass.gtf`
* `rmsk_TEIndividual.gtf`

You can make a gtf with both genes and TEs in teh sam file by using the code below. **However, you will get more TEs counts if you run featureCounts with genes and TEs in separate files**

```bash
# Gene + Individual TEs
cat [gene.annotation.gtf] rmsk_TEIndividual.gtf | sort -k1,1 -k4,4n > TEIndividual_and_genes.gtf

# Gene + Class TEs
cat [gene.annotation.gtf] rmsk_TEClass.gtf | sort -k1,1 -k4,4n > TEClass_and_genes.gtf
```

### STEP 4.1 - Get the counts of the genes and TEs at the individual level

For this we only used the uniquelly mapped bams

```bash
# You can easily get a list of bam files with
bam_files=$(ls ${bam_folder}*.bam)

# For TEs paired-end
featureCounts -a "rmsk_TEIndividual.gtf" \ # Annotation file individual TEs
-o "${output_folder}/Counts_TEIndividual.txt" \
-t exon \
-g gene_name \
-p -C \
-T 16 \
${bam_files} # Example "sample1.bam sample2.bam sample3.bam" 

# For genes paired-end
featureCounts -a "gencode.v44.annotation_gene.gtf" \ # Annotation file for genes
-o "${output_folder}/Counts_Genes_Unique.txt" \
-t gene \
-g gene_name \
-p -C \
-T 16 \
${bam_files} # Example "sample1.bam sample2.bam sample3.bam" 
```

### STEP 4.2 - Get the counts of the genes and TEs at the subfamily level

For this we only used the multi-mapped bams

```bash
# You can easily get a list of bam files with
bam_files=$(ls ${bam_folder}*.bam)

# For TEs paired-end
featureCounts -a "rmsk_TEClass.gtf" \ # Annotation file with families of TEs
-o "${output_folder}/Counts_TEClass.txt" \
-t exon \
-g gene_name \
-M -p -C \
-T 16 \
${bam_files} # Example "sample1.bam sample2.bam sample3.bam" 
```

## Just Press Enter

The script `TE_pipeline.sh` does all of these steps in one go! You just have to give it the fasta_folder, the output folders, the threads you want to use, the mode (unique or multi) and the end-type (single [SE] or paired [PE]).

```bash
bash TE_pipeline.sh [fasta_folder] [output_bams_folder] [output_counts_folder] [threads] [mode] [end]
```

>**Important**
>
> To run `TE_pipeline.sh` you should define some constants that will probably be the same every time you use this script, but the path may change between people. So the first time you run, you should change the paths in lines 23-26.

## STAR Arguments explained

For more check: [STAR Manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

* `--alignIntronMax 500000`: Maximum Intro size
* `--alignMatesGapMax 500000`: Maximum gap size between two mates
* `--alignEndsType EndToEnd`: force end-to-end read alignment, do not soft-clip
* `--outFilterMultiNmax 1`: max number of multiple alignments allowed for a read. If exceeded, the read is considered unmapped
* `--outFilterMismatchNmax 999`: maximum number of mismatched per pair, large number of switches off this filter
* `--outFilterMismatchNoverReadLmax 0.04`: max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.04*200=8 for the paired read
* `--outSAMmultNmax 1`: parameter limits the number of output alignments (SAM lines) for multimappers
* `--seedMultimapNmax`: 20000: only pieces that map fewer than this value are utilized in the stiching procedure

## FeatureCounts Arguments explained

For more check: [featureCounts Docs](https://subread.sourceforge.net/featureCounts.html)

* `-a`: Annotation file
* `-o`: output file
* `-g`: Specify attribute type in a GTF annotation used for counting, default: gene_id
* `-s`: strand. 0 - unstranded; 1- stranded; 2- reversly stranded
* `-C`: don’t count read pairs that have their two ends mapping to different chromosomes or mapping to different strands
* `-p`: paired-end reads. Remove for single-end
* `-M`: Multi-mapping reads will also be counted
