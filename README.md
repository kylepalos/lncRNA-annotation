# lncRNA-annotation

This repo will organize the work being done for the Arabidopsis re-annotation project where Andrew, Caylyn, and I are going to re-map high confidence datasets for lncRNA annotation.

Software versions used:
```
Automatic read quality and adapter trimming: fastp v0.23.2
Read quality detection: FastQC v0.12.1
Read quality collection and reporting: multiQC v1.14
Detect strandedness of library: salmon v1.10.0
Map reads to genome: Hisat2 v2.2.1
Assemble transcripts: StringTie v2.2.1
```

Basic command used for each program and justification:
To start, sequencing experiments will be separated into single-end or paired-end directories so that bash scripting is easier.

Automatic read trimming:
fastp single end:

```
Fastp -i input.fastq.gz -o trimmed_directory/input_trimmed.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```
