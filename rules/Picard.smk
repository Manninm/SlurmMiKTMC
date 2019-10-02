rule Picard:
	input:
		"{sample}/{sample}Aligned.sortedByCoord.out.bam"
	output:
		"{sample}/{sample}picardresults.txt"
	log: "logs/{sample}Picard.log"
	params:
		outprefix = "{sample}/{sample}",

	shell:
		"""
		java -jar {PicardPath} CollectRnaSeqMetrics REF_FLAT={reFlat} STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND CHART_OUTPUT={params.outprefix}coverage.pdf INPUT={input} OUTPUT={output}
		"""
