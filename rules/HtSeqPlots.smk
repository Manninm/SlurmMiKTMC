rule HtSeqBoxPlot:
	input:
		"HtSeqCounts/CountsGt0_voom.txt"
	output:
		"Plots/CountsGt0_voomBoxPlot.png"
	log: "logs/HtSeqBox.log"
	shell:
		"""
		cd HtSeqCounts/; Rscript ../scripts/GeneBoxPlot.R CountsGt0_voom.txt
		"""
rule HtSeqPca:
	input:
		"HtSeqCounts/CountsGt0_voom_filtered.txt"
	output:
		"Plots/CountsGt0_voom_filteredPCA1v2&2v3&3v4.pdf"
	log: "logs/HtSeqPca.log"
	shell:
		"""
		cd HtSeqCounts/; Rscript ../scripts/GenePca.R CountsGt0_voom_filtered.txt
		"""
rule HtSeqCluster:
	input:
		"HtSeqCounts/CountsGt0_voom_filtered.txt"
	output:
		"Plots/CountsGt0_voom_filtered.Ward.D.{GeneBootStraps}.pdf",
		"Plots/CountsGt0_voom_filtered.Ward.D2.{GeneBootStraps}.pdf"
	log: "logs/HtSeqCluster{GeneBootStraps}.log"
	threads: HcCores
	shell:
		"""
		cd HtSeqCounts/; Rscript ../scripts/GeneCluster.R CountsGt0_voom_filtered.txt {threads} {GeneBootStraps} 
		"""