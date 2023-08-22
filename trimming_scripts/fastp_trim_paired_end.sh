for i in *_1.fastq.gz
 do
 name=$(basename ${i} _1.fastq.gz);
 fastp -i ${name}_1.fastq.gz -I ${name}_2.fastq.gz \
 -o $HOME/trim_fastq/paired_end/${name}_1.fastq.gz \
 -O $HOME/trim_fastq/paired_end/${name}_2.fastq.gz \
 -h $HOME/trim_reports/paired_end/${name}.html \
 -j $HOME/trim_reports/paired_end/${name}.json \
 -w 16
 done