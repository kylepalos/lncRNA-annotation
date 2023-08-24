#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 input_file_path"
    exit 1
fi

input_file="$1"

# Loop through each line in the input file
while IFS=$'\t' read -r col1 col2; do
    if [ "$col2" = "SF" ]; then
        basename="$col1"
        name=$(basename "$col1" .fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -U "${name}.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --rna-strandness F \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_single.txt" \
        -p 50 | samtools sort -@ 50 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/se_sf/${basename}.bam"
    elif [ "$col2" = "SR" ]; then
        basename="$col1"
        name=$(basename "$col1" .fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -U "${name}.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --rna-strandness R \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_single.txt" \
        -p 50 | samtools sort -@ 50 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/se_sr/${basename}.bam"
    elif [ "$col2" = "U" ]; then
        basename="$col1"
        name=$(basename "$col1" _1.fastq.gz)
        hisat2 \
        -x /home/kpalos/new_arabidopsis_lncRNAs/genome/t2t_for_mapping/at_t2t \
        -U "${name}.fastq.gz" \
        --max-intronlen 1500 \
        --dta-cufflinks \
        --summary-file "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/summary_files/${basename}_summary_single.txt" \
        -p 50 | samtools sort -@ 50 -o "/home/kpalos/new_arabidopsis_lncRNAs/mapped_reads/se_u/${basename}.bam"
    fi
done < "$input_file"
