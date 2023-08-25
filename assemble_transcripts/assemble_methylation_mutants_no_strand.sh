# for assembling transcripts from the methylation mutant
# Do not specify strand
# uncertain whether StringTie actually benefits from strand specification
for i in *.bam
	do
	name=$(basename ${i} .bam);
	stringtie ${name}.bam \
	-G athal_t2t_liftover_gffread.gtf \
	-o $HOME/transcript_assembly/methylation_mutants/specify_strand/${name}.gtf \
	-f 0.05 -j 5 -c 5 -s 15 -p 40;
	done