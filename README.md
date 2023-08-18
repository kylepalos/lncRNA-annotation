# lncRNA-annotation

This repo will organize the work being done for the Arabidopsis re-annotation project where Andrew, Caylyn, and I are going to re-map high confidence datasets for lncRNA annotation.

**Software versions used:**
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

**Automatic read trimming:**

fastp single-end:

```
Fastp -i input.fastq.gz -o trimmed_directory/input_trimmed.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```

fastp paired -end:

```
Fastp -i input_1.fastq.gz -I input_2.fastq.gz \
-o trimmed_directory/input_trimmed_1.fastq.gz -O input_trimmed_2.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```

Automatic adapter and quality detection followed by trimming will avoid manual and error-prone trimming (phrase better)



**Library strandedness detection:**

```
# make the index from Araport11 transcripts
salmon index -t at_transcripts.fa -i athal_index

# quantify the reads
# overkill but we'll scrape the library type from the log file
# single end
salmon quant -l A -r input_trimmed.fastq.gz \
-i athal_index \
-o salmon_detection/input_quant

# paired end
salmon quant -l A -1 input_trimmed_1.fastq.gz -2 input_trimmed_2.fastq.gz \
-i athal_index \
-o salmon_detection/input_quant
```

Salmon is one of the only programs I know that so easily reports the strandedness of a dataset.

Include the script to scrape the log file for library type detection.
