for i in *_1.fastq.gz
 do
 name=$(basename ${i} _1.fastq.gz); # establish basename
 fastp -i ${name}_1.fastq.gz -I ${name}_2.fastq.gz \ # input forward and reverse mate
 -o /$HOME/trim_fastq/methylation_mutants/${name}_1.fastq.gz \ # output forward mate
 -O /$HOME/trim_fastq/methylation_mutants/${name}_2.fastq.gz \ # output reverse mate
 -h /$HOME/trim_reports/methylation_mutants/${name}.html \ # output html report
 -j /$HOME/trim_reports/methylation_mutants/${name}.json \ # output json report
 -w 16 # 16 threads
 done