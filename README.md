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
Classify transcripts: gffcompare v0.12.6
```

**The remainder of this readme is currently an overview of the commands used for each program and some justification:**
To start, sequencing experiments will be separated into single-end or paired-end directories so that bash scripting is easier.

Additionally, experiments dealing with DNA methylation mutant experiments will be kept sepate for the entire workflow.

**Automatic read trimming:**

fastp single-end:

```
fastp -i input.fastq.gz -o trimmed_directory/input_trimmed.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```

fastp paired -end:

```
fastp -i input_1.fastq.gz -I input_2.fastq.gz \
-o trimmed_directory/input_trimmed_1.fastq.gz -O input_trimmed_2.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```

Automatic adapter and quality detection followed by trimming will avoid manual and error-prone processing using other tools where we have to manually enter trimming parameters.



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

Caylyn pushed this script to the repo, "extract_library_type.sh"


**Read mapping prep:**

```
# Make the Hisat2 index from Araport11 GTF
# Obviously change with new version for T2T genome
# This will improve anchoring during read mapping

extract_exons.py Araport11.gtf > Athaliana.exon

extract_splice_sites.py Araport11.gtf > Athaliana.ss

# Build command:

# will put new Athaliana genome fasta version here
hisat2-build --ss Athaliana.ss --exon Athaliana.exon Arabidopsis_thaliana.TAIR10.dna.toplevel.fa AT_hisat -p 32

```

**Read mapping:**

```
hisat2 -x AT_hisat \ #index

--max-intronlen 1500 \ # <1% of introns are greater than 1kb, see: 10.3390/genes8080200

--dta-cufflinks \ # recommended for conservative transcript assembly, see https://github.com/gpertea/stringtie/issues/204#issuecomment-435409259

--rna-strandness ${} \ # as reported by Salmon and converted to the proper Hisat2 argument

-U $trimmed_reads.fastq.gz \ # for single end

-1 $trimmed_reads_forward_1.fastq.gz -2 $trimmed_reads_reverse_2.fastq.gz | # for paired end, obviously only use single or paired end

samtools sort -o $file_name.bam - # pipe output to bam conversion
```


**Transcript assembly:**

```
# first index all the bam files
ls *.bam | parallel samtools index '{}'

# assemble methylation mutant samples separately from all other samples
# assemble paired and single end samples to the same output directory because they will be merged together

stringtie ${name}.bam \

-G Araport_reference.gtf \ # specify the reference annotation to guide new assemblies

-o ${name}_stringtie.gtf \ # output

-f 0.05 \ # alternative isoforms must be present at > 5% abundance of the primary isoform
# Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts.

-j 5 \ # the number of spliced reads that should align across a junction (junction coverage) - reduces spurious spliced transcripts

-c 5 \ # minimum read coverage allowed for the predicted transcripts

-s 15  \# minimum read coverage allowed for single-exon transcripts 

[--rf, --fr, or omit if unstranded]
```

**Merge assemblies:**

```
stringtie --merge [all single and paired-end outputs].gtf \
-G reference_annotation.gtf \
-o merged_output.gtf
```




**Transcript classification:**

```
gffcompare -r reference_annotation.gtf -o Athaliana_gffcompare merged_output.gtf
```

Gffcompare generates class codes for all newly assembled transcripts, [see here](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

Class code 'u' are new intergenic transcripts (may eventually become lincRNAs)

Class code 'x' are antisense transcripts to annotated genes (may eventually become antisense lncRNAs)

**Transcripts from these classes will be fed into:**
1. [CPC2](http://cpc2.gao-lab.org/run_cpc2_program.php) - coding classification
2. [rFam](https://rfam.org/search) - similarity to "housekeeping" RNAs (snoRNA, rRNA, tRNA, etc.)
3. [PfamScan](https://www.ebi.ac.uk/Tools/pfa/pfamscan/) - similarity to known protein domains


**Liftover of Araport11 annotation to new T2T genome:**

**Note: The soon to be published annotation will actually be more manually and accurately curated re-annotation for the T2T genome. I am doing this general liftover for antisense lncRNA identification purposes onyl.**

I will be using the [Liftoff](https://github.com/agshumate/Liftoff) tool coming out of Steven Salzberg's group, which uses Minimap2 to accurately map annotations between assemblies.

I am using a slightly custom Araport11 gff3 annotation that has been cleaned using [AGAT's](https://github.com/NBISweden/AGAT) agat_convert_sp_gxf2gxf.pl script for annotation standardization.

Additionally, I will provide a list of feature types beyond what Liftoff natively recognizes, this is essentially just a .txt file that says:

```
lnc_RNA
miRNA
tRNA
ncRNA
snoRNA
snRNA
rRNA
```

provided through the -f argument.

The command for liftover was:
```
python3 liftoff/build/lib/liftoff/run_liftoff.py \
athal_t2t.fa tair_no_organelle.fa \ # target then reference
-g athaliana_araport_no_organelle.gff3 \ # annotation of features to lift over
-f feature_types_file.txt \ # the feature file above
-o athal_t2t_liftover.gff3 # the output annotation
```
