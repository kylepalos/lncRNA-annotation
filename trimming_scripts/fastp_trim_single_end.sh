for i in *.fastq.gz
 do
 name=$(basename ${i} .fastq.gz); # isolate the basename
 fastp -i ${name}.fastq.gz -o $HOME/trim_fastq/single_end/${name}.fastq.gz \ # input fastq and output trimmed fastq
 -h $HOME/trim_reports/single_end/${name}.html \ # where to output the HTML reports
 -j $HOME/trim_reports/single_end/${name}.json \ # where to output the JSON reports
 -w 8 # number of threads
 done