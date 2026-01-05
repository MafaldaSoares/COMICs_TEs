#!/bin/bash
#########################################################
# Script to align and count the reads for the TE pipeline
# Date: 29-10-2025
#########################################################

set -euo pipefail

################# USAGE AND ARGUMENTS ###################

if [[ $# -lt 6 ]]; then
  echo "Usage: $0 <fasta_folder> <output_bam_folder> <output_counts_folder> <threads> <mode: unique|multi> <end: SE|PE>"
  echo "Example: $0 Fasta_Folder/ BAM_Folder/ Counts_Folder/ 32 unique PE"
  exit 1
fi

samples_folder="$1"
bam_folder="$2"
counts_folder="$3"
threads="$4"
mode="$5"         # "unique" or "multi"
end="$6"       # "SE" (single-end) or "PE" (paired-end)
genomeDir="GRCh38.primary_assembly/STAR_index/" # STAR_Index folder path
TEIndividual_gtf="rmsk_TEIndividual.gtf" # rmsk_TEIndividual.gtf path
TEClass_gtf="rmsk_TEClass.gtf" #rmsk TEClass.gtf path
Gene_gtf="GRCh38.primary_assembly/GeneAnnotations/gencode.v44.annotation_gene.gtf" #gene.annotation.gtf path


######################### SETUP #########################

mkdir -p "${bam_folder}"
ulimit -n 10000

echo "-------------------------------------------"
echo "Input folder:   ${samples_folder}"
echo "Output folder:  ${bam_folder}"
echo "Threads:        ${threads}"
echo "Mode:           ${mode}"
echo "End:            ${end}"
echo "Genome index:   ${genomeDir}"
echo "-------------------------------------------"

# Unique vs Multimapping settings
if [[ "${mode}" == "unique" ]]; then
  outFilterMultimapNmax=1
  outFilterMismatchNoverLmax=0.03
  extra_align_args="--alignIntronMax 500000 --alignMatesGapMax 500000 --alignEndsType EndToEnd"
else
  outFilterMultimapNmax=5000
  outFilterMismatchNoverLmax=0.06
  extra_align_args=""
fi

##################### STAR ALIGNMENT ####################

echo "Running STAR alignment..."

# For paired end cases
if [[ "${end}" == "PE" ]]; then
  files=(${samples_folder}*_1.fq*)
  total_samples=${#files[@]}
  counter=1

  for file in "${files[@]}"; do
    filename=$(basename "${file}")
    sample="${filename%_*}"
    echo "(${counter}/${total_samples}) Mapping paired-end: ${sample}"

    STAR --runThreadN "${threads}" \
      --genomeDir "${genomeDir}" \
      --limitSjdbInsertNsj 6000000 \
      ${extra_align_args} \
      --outFileNamePrefix "${bam_folder}/${sample}_" \
      --outFilterMultimapNmax "${outFilterMultimapNmax}" \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax "${outFilterMismatchNoverLmax}" \
      --seedMultimapNmax 20000 \
      --outSAMattributes All \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesIn "${file}" "${file%_*}_2.fq" 

    ((counter++))
  done


# For single end
elif [[ "${end}" == "SE" ]]; then
  files=(${samples_folder}*.fq*)
  total_samples=${#files[@]}
  counter=1

  for file in "${files[@]}"; do
    filename=$(basename "${file}")
    sample="${filename%.*}"
    echo "(${counter}/${total_samples}) Mapping single-end: ${sample}"

    STAR --runThreadN "${threads}" \
      --genomeDir "${genomeDir}" \
      ${extra_align_args} \
      --outFileNamePrefix "${bam_folder}/${sample}_" \
      --outFilterMultimapNmax "${outFilterMultimapNmax}" \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax "${outFilterMismatchNoverLmax}" \
      --seedMultimapNmax 20000 \
      --outSAMattributes All \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesIn "${file}" 

    ((counter++))
  done
else
  echo "Error: End must be 'PE' (paired-end) or 'SE' (single-end)"
  exit 1
fi

################## STAR SUMMARY EXTRACTION #################

echo "Extracting STAR alignment summaries..."

summary_file="${bam_folder}/STAR_${mode}_alignment_summary.tsv"
echo -e "Sample\tInput_Reads\tUniquely_Reads\tUniquely_Percent\tMulti_Reads\tMulti_Percent\tUnmapped_TooManyMismatch\tUnmapped_TooManyMismatch_Percent\tUnmapped_TooShort\tUnmapped_TooShort_Percent\tUnmapped_Other\tUnmapped_Other_Percent" > "$summary_file"

for logfile in "${bam_folder}"/*_Log.final.out; do
  [[ -f "$logfile" ]] || continue
  sample=$(basename "$logfile" | sed 's/_Log\.final\.out$//')

  input_reads=$(grep "Number of input reads" "$logfile" | awk '{print $NF}')
  uniquely_mapped_reads=$(grep "Uniquely mapped reads number" "$logfile" | awk '{print $NF}')
  uniquely_mapped_percent=$(grep "Uniquely mapped reads %" "$logfile" | awk '{print $NF}')
  multi_mapped_reads=$(grep "Number of reads mapped to multiple loci" "$logfile" | awk '{print $NF}')
  multi_mapped_percent=$(grep "% of reads mapped to multiple loci" "$logfile" | awk '{print $NF}')
  unmapped_mismatch=$(grep "Number of reads unmapped: too many mismatches" "$logfile" | awk '{print $NF}')
  unmapped_mismatch_percent=$(grep "% of reads unmapped: too many mismatches" "$logfile" | awk '{print $NF}')
  unmapped_short=$(grep "Number of reads unmapped: too short" "$logfile" | awk '{print $NF}')
  unmapped_short_percent=$(grep "% of reads unmapped: too short" "$logfile" | awk '{print $NF}')
  unmapped_other=$(grep "Number of reads unmapped: other" "$logfile" | awk '{print $NF}')
  unmapped_other_percent=$(grep "% of reads unmapped: other" "$logfile" | awk '{print $NF}')

  echo -e "${sample}\t${input_reads}\t${uniquely_mapped_reads}\t${uniquely_mapped_percent}\t${multi_mapped_reads}\t${multi_mapped_percent}\t${unmapped_mismatch}\t${unmapped_mismatch_percent}\t${unmapped_short}\t${unmapped_short_percent}\t${unmapped_other}\t${unmapped_other_percent}" >> "$summary_file"
done

echo "Summary written to: ${summary_file}"
echo "All done!"


#################### FEATURE COUNTS ####################

echo "Running featureCounts quantification..."

# Collect all BAM files generated by STAR
bam_files=$(ls ${bam_folder}*.bam 2>/dev/null || true) # Find bam files and don't crash if none are found

if [[ -z "$bam_files" ]]; then
  echo "ERROR: No BAM files found in ${bam_folder}. Skipping featureCounts."
  exit 1
fi

# Adjust the number of threads for featurecounts, since it can only take 64 maximum
if (( threads > 64 )); then
  fc_threads=32
  echo "featureCounts threads was limited to 32 threads for safety (input was ${threads})"
else
  fc_threads="${threads}"
fi

# Select annotation and featureCounts parameters based on mode/end
if [[ "${mode}" == "unique" && "${end}" == "PE" ]]; then
  # Arguments for TE
  annotation="${TEIndividual_gtf}"
  outfile="${counts_folder}/Counts_TE_Individual_AllSamples.txt"
  fc_args="-t exon -g gene_name -p -C"
  # Arguments for gene 
  annotation_gene="${Gene_gtf}"
  outfile_gene="${counts_folder}/Counts_Gene_Individual_AllSamples.txt"
  fc_args_gene="-t gene -g gene_name -p -C"

elif [[ "${mode}" == "unique" && "${end}" == "SE" ]]; then
  annotation="${TEIndividual_gtf}"
  outfile="${counts_folder}/Counts_TE_Individual_AllSamples.txt"
  fc_args="-t exon -g gene_name"
  # Arguments for gene 
  annotation_gene="${Gene_gtf}"
  outfile_gene="${counts_folder}/Counts_Gene_Individual_AllSamples.txt"
  fc_args_gene="-t gene -g gene_name"

elif [[ "${mode}" == "multi" && "${end}" == "PE" ]]; then
  annotation="${TEClass_gtf}"
  outfile="${counts_folder}/Counts_TE_Class_AllSamples.txt"
  fc_args="-t exon -g gene_name -M -p -C"

elif [[ "${mode}" == "multi" && "${end}" == "SE" ]]; then
  annotation="${TEClass_gtf}"
  outfile="${counts_folder}/Counts_TE_Class_AllSamples.txt"
  fc_args="-t exon -g gene_name -M"

else
  echo "ERROR: Invalid combination of mode (${mode}) and end (${end})."
  exit 1
fi


# Run featureCounts
if [[ "${mode}" == "unique" ]]; then
  
  # Count TEs
  echo "-------------------------------------------"
  echo "featureCounts TE parameters:"
  echo "  Annotation: ${annotation}"
  echo "  Output:     ${outfile}"
  echo "  Threads:    ${fc_threads}"
  echo "  Extra args: ${fc_args}"
  echo "-------------------------------------------"
  
  featureCounts -a "${annotation}" \
    -o "${outfile}" \
    ${fc_args} \
    -T "${fc_threads}" \
    ${bam_files}
  
  # Count Genes
  echo "-------------------------------------------"
  echo "featureCounts Gene parameters:"
  echo "  Annotation: ${annotation_gene}"
  echo "  Output:     ${outfile_gene}"
  echo "  Threads:    ${fc_threads}"
  echo "  Extra args: ${fc_args}_gene"
  echo "-------------------------------------------"
  
  featureCounts -a "${annotation_gene}" \
    -o "${outfile_gene}" \
    ${fc_args_gene} \
    -T "${fc_threads}" \
    ${bam_files}

elif [[ "${mode}" == "multi" ]]; then
  
  # Count TEs
  echo "-------------------------------------------"
  echo "featureCounts TE parameters:"
  echo "  Annotation: ${annotation}"
  echo "  Output:     ${outfile}"
  echo "  Threads:    ${fc_threads}"
  echo "  Extra args: ${fc_args}"
  echo "-------------------------------------------"
  
  featureCounts -a "${annotation}" \
    -o "${outfile}" \
    ${fc_args} \
    -T "${fc_threads}" \
    ${bam_files}
fi

echo "featureCounts completed successfully!"
echo "Counts file written to: ${outfile}"
