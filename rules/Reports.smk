rule StarCollect:
    input:
        input=expand('{sample}/{sample}Log.final.out',sample=SAMPLES)
    output:	
        'Reports/averageMappedLength.txt',
        'Reports/input.txt',
        'Reports/mappedPercent.txt',
        'Reports/percentMultimappers.txt',
        'Reports/unmappedTooShort.txt'
    shell:
        """
		for i in $( find . -maxdepth 2 -name "*Log.final.out" ); do echo $i >> Reports/input.txt; cat $i | grep "Number of input reads" >> Reports/input.txt; done;
		for i in $( find . -maxdepth 2 -name "*Log.final.out" ); do echo $i >> Reports/averageMappedLength.txt; cat $i | grep "Average mapped length" >> Reports/averageMappedLength.txt; done;
		for i in $( find . -maxdepth 2 -name "*Log.final.out" ); do echo $i >> Reports/percentMultimappers.txt; cat $i | grep "% of reads
		mapped to multiple loci" >> Reports/percentMultimappers.txt; done;
		for i in $( find . -maxdepth 2 -name "*Log.final.out" ); do echo $i >> Reports/unmappedTooShort.txt; cat $i | grep " % of reads
		unmapped: too short" >> Reports/unmappedTooShort.txt; done
		for i in $( find . -maxdepth 2 -name "*Log.final.out" ); do echo $i >> Reports/mappedPercent.txt; cat $i | grep "Uniquely mapped
		reads %" >> Reports/mappedPercent.txt; done
		"""
rule PicardGather:
	input:
		expand("{sample}/{sample}picardresults.txt",sample=SAMPLES)
	output:
		"Reports/picard_classifications.txt"
	shell:
		"""
		for i in $( find . -maxdepth 2 -wholename "*picardresults.txt" -type f ); do echo $i >> Reports/picard_classifications.txt; head -n 8 $i |tail -n 2 >> Reports/picard_classifications.txt; done
		"""
rule MultiQc:
	input:
		expand("{sample}/{sample}picardresults.txt",sample=SAMPLES),
		expand("{sample}/{sample}Log.final.out",sample=SAMPLES),
		expand("{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/summary.txt",sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/summary.txt",sample=SAMPLES),
		expand("{sample}/{sample}_htseq.cnt",sample=SAMPLES),
		expand("{sample}/{sample}_screen.txt",sample=SAMPLES)
	output:
		"multiqc_report.html",
	shell:
		"""
		multiqc --module htseq --module star --module picard --module fastqc --module fastq_screen .
		"""		