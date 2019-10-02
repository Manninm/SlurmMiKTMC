rule String_Tie:
	input:
		"{sample}/{sample}Aligned.sortedByCoord.out.bam"
	output:
		"{sample}/{sample}_onlyKnown.gtf"
	log: "logs/{sample}StringTie.log"
	threads: 4
	shell:
		"""
		stringtie -p {threads} -eB -G {MYGTF} -o {output} {input}
		"""
