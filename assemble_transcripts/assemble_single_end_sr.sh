# # for assembling transcripts from single end samples
# with SR library type reported by Salmon
for i in *.bam
	do
	name=$(basename ${i} .bam);
	stringtie ${name}.bam \
	-G athal_t2t_liftover_gffread.gtf \
	-o $HOME/transcript_assembly/single_end/sr/${name}.gtf \
	-f 0.05 -j 5 -c 5 -s 15 --rf -p 40;
done