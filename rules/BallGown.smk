if not config["GroupdFile"]:
	os.system('Rscript scripts/Table.R')
	logger.info("No Group File provided!")
rule Ball_Gown:
	input:
		expand("{sample}/{sample}_onlyKnown.gtf",sample=SAMPLES)
	output:
		"BallGown/gene_expression_table.txt",
		"BallGown/transcript_fpkm.txt",
		"BallGown/transcript_cov.txt",
		"BallGown/whole_tx_table.txt",
		"BallGown/filt_gene_expression_table.txt",
		"BallGown/filt_transcript_fpkm.txt",
		"BallGown/filt_transcript_cov.txt",
		"BallGown/filt_whole_tx_table.txt",
	log: "logs/BallGown.log"
	shell:
		"""
		Rscript scripts/BallGown.R
		"""