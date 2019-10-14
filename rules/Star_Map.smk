rule STAR_fq:
    input: 
        R1 = "{sample}/{sample}.R1.fq.gz",
        R2 = "{sample}/{sample}.R2.fq.gz",
    output: 
        "{sample}/{sample}Aligned.sortedByCoord.out.bam",
        '{sample}/{sample}Log.final.out'
    params:
        jobname = "{sample}",
        outprefix = "{sample}/{sample}"
    threads: StarThreads
	message: "aligning {input} using STAR: {threads} threads"
	log: "logs/{sample}.starlog"
	shell:
		"""
        STAR --twopassMode Basic --runThreadN {threads} \
        --genomeDir {STARINDEX} --outSAMtype BAM SortedByCoordinate \
        --outBAMcompression 10 --outSAMstrandField intronMotif \
        --outBAMsortingThreadN {SortThreads} \
        --bamRemoveDuplicatesType UniqueIdentical --readFilesCommand gunzip -c --readFilesIn {input.R1} {input.R2} \
        --outFileNamePrefix {params.outprefix}
		"""


