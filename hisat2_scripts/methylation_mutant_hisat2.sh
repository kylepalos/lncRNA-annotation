#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 input_file_path"
    exit 1
fi

input_file="$1"

# Loop through each line in the input file
while IFS=$'\t' read -r col1 col2; do
    if [ "$col2" = "ISF" ]; then
        basename="$col1"
        name=$(basename "$col1" _1.fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -1 "${name}_1.fastq.gz" \
        -2 "${name}_2.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --rna-strandness FR \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_methylation_mutant.txt" \
        -p 80 | samtools sort -@ 80 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/methylation_mutants/${basename}.bam"
    elif [ "$col2" = "ISR" ]; then
        basename="$col1"
        name=$(basename "$col1" _1.fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -1 "${name}_1.fastq.gz" \
        -2 "${name}_2.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --rna-strandness RF \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_methylation_mutant.txt" \
        -p 80 | samtools sort -@ 80 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/methylation_mutants/${basename}.bam"
    elif [ "$col2" = "IU" ]; then
        basename="$col1"
        name=$(basename "$col1" _1.fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -1 "${name}_1.fastq.gz" \
        -2 "${name}_2.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_methylation_mutant.txt" \
        -p 80 | samtools sort -@ 80 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/methylation_mutants/${basename}.bam"
    fi
done < "$input_file"
