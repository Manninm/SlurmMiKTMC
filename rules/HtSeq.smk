rule HtSeq_Count:
	input: "{sample}/{sample}Aligned.sortedByCoord.out.bam"
	output: "{sample}/{sample}_htseq.cnt"
	log: "logs/{sample}HtSeq.log"
	params: 
		jobname = "{sample}"
	threads: 2
	message: "htseq-count {input} : {threads} threads"
	shell:
		"""
		htseq-count -f bam --stranded=no --idattr=gene_id --order=pos -m union {input} {MYGTF} > {output}
		"""
rule HtSeq_Move:
	input:
		"{sample}/{sample}_htseq.cnt"
	output:
		"HtSeqCounts/{sample}_htseq.cnt"
	shell:
		"""
		scripts/HtSeqRename.sh
		"""
rule HtSeq_Matrix:
	input:
		expand("HtSeqCounts/{sample}_htseq.cnt",sample=SAMPLES)
	output:
		"HtSeqCounts/allCounts.txt",
		"HtSeqCounts/CountsGt0.txt",
		"HtSeqCounts/CountsGt0_voom.txt",
		"HtSeqCounts/CountsGt0_voom_filtered.txt"
	shell:
		"cd HtSeqCounts/; Rscript ../scripts/HtSeqVoomMerge.R"