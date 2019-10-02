rule BallGownGeneBoxPlot:
	input:
		"BallGown/filt_gene_expression_table.txt"
	output:
		"Plots/filt_gene_expression_tableBoxPlot.png"
	log: "logs/GeneBox.log"
	shell:
		"""
		cd BallGown/; Rscript ../scripts/GeneBoxPlot.R filt_gene_expression_table.txt
		"""
rule BallGownGenePca:
	input:
		"BallGown/filt_gene_expression_table.txt"
	output:
		"Plots/filt_gene_expression_tablePCA1v2&2v3&3v4.pdf"
	log: "logs/GenePca.log"
	shell:
		"""
		cd BallGown/; Rscript ../scripts/GenePca.R filt_gene_expression_table.txt
		"""
rule BallGownGeneCluster:
	input:
		"BallGown/filt_gene_expression_table.txt"
	output:
		"Plots/filt_gene_expression_table.Ward.D.{GeneBootStraps}.pdf",
		"Plots/filt_gene_expression_table.Ward.D2.{GeneBootStraps}.pdf"
	log: "logs/GeneCluster{GeneBootStraps}.log"
	threads: HcCores
	shell:
		"""
		cd BallGown/; Rscript ../scripts/GeneCluster.R filt_gene_expression_table.txt {threads} {GeneBootStraps} 
		"""