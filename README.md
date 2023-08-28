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

Additionally, experiments dealing with DNA methylation mutant experiments will be processed but kept separate for the majority of the workflow.


## Automatic read trimming:

Automatic adapter and quality detection followed by trimming will avoid manual and error-prone processing using other tools where we have to manually enter trimming parameters.


**fastp single-end:**

```
fastp -i input.fastq.gz -o trimmed_directory/input_trimmed.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```

**fastp paired-end:**

```
fastp -i input_1.fastq.gz -I input_2.fastq.gz \
-o trimmed_directory/input_trimmed_1.fastq.gz -O input_trimmed_2.fastq.gz \
-h html_reports/input.html \
-j json_reports/input.json
```




## Library strandedness detection:

Salmon is one of the only programs I know that so easily reports the strandedness of a dataset.

We'll make an output file that maps every SRA ID to a strandedness code for proper read mapping and transript assembly later on.

Assess strandedness of data using Araport11 transcripts.

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


## Liftover for mapping:

**Note: The soon to be published annotation will actually be a more manually and accurately curated re-annotation for the T2T genome. I am doing this general liftover for lncRNA identification purposes only.**

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



## Read mapping prep:

We're going to use Hisat2 because StringTie was developed with Hisat2 formatted alignment files in mind, see [this StringTie GitHub issue for more context](https://github.com/gpertea/stringtie/issues/204).

**Making the Hisat2 index for mapping:**

```
# first convert gff to gtf so Hisat can make relevant files
gffread -T athal_t2t_liftover.gff3 > athal_t2t_liftover_gffread.gtf

# then use Hisat2 extract exons and splice sites so we can make the index:

python3 extract_exons.py athal_t2t_liftover_gffread.gtf > at_t2t.exon

python3 extract_splice_sites.py athal_t2t_liftover_gffread.gtf > at_t2t.ss

hisat2-build --ss at_t2t.ss --exon at_t2t.exon athal_t2t.fa at_t2t -p 128
```


## Read mapping:

Going to use relatively simple parameters, each argument is commented for context.

```
hisat2 -x AT_hisat \ #index

--max-intronlen 1500 \ # <1% of introns are greater than 1kb, see: 10.3390/genes8080200

--dta-cufflinks \ # recommended for conservative transcript assembly, see https://github.com/gpertea/stringtie/issues/204#issuecomment-435409259

--rna-strandness ${} \ # as reported by Salmon and converted to the proper Hisat2 argument

-U $trimmed_reads.fastq.gz \ # for single end

-1 $trimmed_reads_forward_1.fastq.gz -2 $trimmed_reads_reverse_2.fastq.gz | # for paired end, obviously only use single or paired end

samtools sort -o $file_name.bam - # pipe output to bam conversion
```


## Transcript assembly:

Using StringTie to assemble new transcripts from alignment files. Will use liftover annotation as a reference for assembly.

```
# first index all the bam files
ls *.bam | parallel samtools index '{}'

# Run StringTie

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


## Merge assemblies:

Potentially contentious, but I am assembling the methylation mutant samples along with the "regular" samples. This will allow me to quantify the RNA-seq reads against all these new transcripts and I can do post-processing to find transcripts only expressed in these mutant samples. They can then be removed/set aside.

```
stringtie --merge [all single and paired-end outputs].gtf \
-G reference_annotation.gtf \
-o merged_output.gtf
```

**Remove any transcripts which do not get a '+' or '-' designation for strand**

**These are single-exon transcripts which I find are typically unreliable and problematic to have in a gff**

```
awk '$7 != "."' merged_output.gtf > merged_output_no_unstranded.gtf
```



## Transcript classification:

```
gffcompare -r reference_annotation.gtf -o Athaliana_gffcompare merged_output_no_unstranded.gtf
```

Gffcompare generates class codes for all newly assembled transcripts, [see here](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

Class code 'u' are new intergenic transcripts (may eventually become lincRNAs)

Class code 'x' are antisense transcripts to annotated genes (may eventually become antisense lncRNAs)

Class code 'i' are fully intronic transcripts to annotated genes. I will extract out those that are on the opposite strand as the reference.


**Transcripts from these classes will be fed into:**
1. [CPC2](http://cpc2.gao-lab.org/run_cpc2_program.php) - coding classification
2. [rFam](https://rfam.org/search) - similarity to "housekeeping" RNAs (snoRNA, rRNA, tRNA, etc.)
3. [PfamScan](https://www.ebi.ac.uk/Tools/pfa/pfamscan/) - similarity to known protein domains


## Expression re-analysis:

After assembly, it will be useful to see how these "new" genes are expressed across the sampled experiments we used.

I did this doing the following steps:

1. Merge the "liftover" annotation with the new assembly of interest (I extracted out the genes with the class codes of interest from the merged annotation).
2. I then made a decoy aware Salmon index, following [this workflow](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
3. And quantified all the trimmed fastq files against this index:

```
salmon quant -l A -i athal_decoy_aware_index \
-1 reads_1.fastq.gz -2 reads_2.fastq.gz \
--gcBias --posBias --seqBias \
-o name_quant
```






