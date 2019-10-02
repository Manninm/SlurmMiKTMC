rule BallGownTxBoxPlot:
	input:
		"BallGown/transcript_fpkm.txt"
	output:
		"Plots/transcript_fpkmBoxPlot.png"
	log: "logs/TxBox.log"
	shell:
		"""
		cd BallGown/; Rscript ../scripts/TxBoxPlot.R transcript_fpkm.txt
		"""
rule BallGownTxPca:
	input:
		"BallGown/transcript_fpkm.txt"
	output:
		"Plots/transcript_fpkmPCA1v2&2v3&3v4.pdf"
	log: "logs/TxPca.log"
	shell:
		"""
		cd BallGown/; Rscript ../scripts/TxPca.R transcript_fpkm.txt
		"""
rule BallGownTxCluster:
	input:
		"BallGown/transcript_fpkm.txt"
	output:
		"Plots/transcript_fpkm.Ward.D2.{TxBootStraps}.pdf",
		"Plots/transcript_fpkm.Ward.D.{TxBootStraps}.pdf"
	log: "logs/TxCluster{TxBootStraps}.log"	
	threads: HcCores
	shell:
		"""
		cd BallGown/; Rscript ../scripts/TxCluster.R transcript_fpkm.txt {threads} {TxBootStraps} 
		"""