rule fastqc:
	input:
		"{sample}/{sample}.R1.fq.gz",
		"{sample}/{sample}.R2.fq.gz"
	output:
		"{sample}/{sample}.R1_fastqc/Images/duplication_levels.png",
		"{sample}/{sample}.R1_fastqc/Images/per_sequence_gc_content.png",
		"{sample}/{sample}.R1_fastqc/Images/sequence_length_distribution.png",
		"{sample}/{sample}.R1_fastqc/Images/per_base_n_content.png",
		"{sample}/{sample}.R1_fastqc/Images/per_tile_quality.png",
		"{sample}/{sample}.R1_fastqc/Images/kmer_profiles.png",
		"{sample}/{sample}.R1_fastqc/Images/per_base_quality.png",
		"{sample}/{sample}.R1_fastqc/Images/per_base_sequence_content.png",
		"{sample}/{sample}.R1_fastqc/Images/adapter_content.png",
		"{sample}/{sample}.R1_fastqc/Images/per_sequence_quality.png",
		"{sample}/{sample}.R2_fastqc/Images/duplication_levels.png",
		"{sample}/{sample}.R2_fastqc/Images/per_sequence_gc_content.png",
		"{sample}/{sample}.R2_fastqc/Images/sequence_length_distribution.png",
		"{sample}/{sample}.R2_fastqc/Images/per_base_n_content.png",
		"{sample}/{sample}.R2_fastqc/Images/per_tile_quality.png",
		"{sample}/{sample}.R2_fastqc/Images/kmer_profiles.png",
		"{sample}/{sample}.R2_fastqc/Images/per_base_quality.png",
		"{sample}/{sample}.R2_fastqc/Images/per_base_sequence_content.png",
		"{sample}/{sample}.R2_fastqc/Images/adapter_content.png",
		"{sample}/{sample}.R2_fastqc/Images/per_sequence_quality.png",
		"{sample}/{sample}.R1_fastqc/summary.txt",
		"{sample}/{sample}.R2_fastqc/summary.txt",
	log: "logs/{sample}FastQc.log"
	threads: 2
	shell:
		"""
		fastqc -t {threads} --extract {input}
		"""
rule fastqScreen:
	input:
		Fast1="{sample}/{sample}.R1.fq.gz",
		Fast2="{sample}/{sample}.R2.fq.gz"
	output:
		output1="{sample}/{sample}.fq.gz",
		output2="{sample}/{sample}_screen.png",
		output3="{sample}/{sample}_screen.txt"
	log: "logs/{sample}FastScreen.log"
	params: 
		outprefix = "{sample}"
	threads: 4
	priority: 3
	shell:
		"""
		cat {input.Fast1} {input.Fast2} > {output.output1} && /home/manninm/Programs/fastq_screen_v0.14.0/fastq_screen --aligner bowtie2 --quiet --force --threads {threads} {output.output1}
		"""
rule PNGRename:
	input:
		expand("{sample}/{sample}.R1_fastqc/Images/duplication_levels.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_sequence_gc_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/sequence_length_distribution.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_base_n_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_tile_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/kmer_profiles.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_base_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_base_sequence_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/adapter_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R1_fastqc/Images/per_sequence_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/duplication_levels.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_sequence_gc_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/sequence_length_distribution.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_base_n_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_tile_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/kmer_profiles.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_base_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_base_sequence_content.png",sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/adapter_content.png", sample=SAMPLES),
		expand("{sample}/{sample}.R2_fastqc/Images/per_sequence_quality.png", sample=SAMPLES),
		expand("{sample}/{sample}_screen.png",sample=SAMPLES)
	output:
		output1="QcImages/adapter_content_pngs.pdf",
		output2="QcImages/duplication_levels_pngs.pdf",
		output3="QcImages/kmer_profiles_pngs.pdf",
		output4="QcImages/per_base_quality_pngs.pdf",
		output5="QcImages/per_base_sequence_content_pngs.pdf",
		output6="QcImages/per_sequence_gc_content_pngs.pdf",	
		output7='Screen/screen_pngs.pdf',
	run:
		if not os.path.exists(MainDir+'QcImages'):
			os.mkdir(MainDir+'QcImages')
		for i in range(len(SAMPLES)):
			R1paths=os.listdir(SAMPLES[i]+'/'+SAMPLES[i]+'.R1_fastqc/Images/')
			for image in range(len(R1paths)):
				print(R1paths[image])
				shutil.copyfile(SAMPLES[i]+'/'+SAMPLES[i]+'.R1_fastqc/Images/'+R1paths[image],'QcImages/'+SAMPLES[i]+'.R1.'+R1paths[image])
		for i2 in range(len(SAMPLES)):
			R2paths=os.listdir(SAMPLES[i2]+'/'+SAMPLES[i2]+'.R2_fastqc/Images/')
			for image2 in range(len(R2paths)):
				print(R2paths[image2])
				shutil.copyfile(SAMPLES[i2]+'/'+SAMPLES[i2]+'.R2_fastqc/Images/'+R2paths[image2],'QcImages/'+SAMPLES[i2]+'.R2.'+R2paths[image2])
		if not os.path.exists(MainDir+'Screen/'):
			os.mkdir(MainDir+'Screen/')				
		for i in range(len(SAMPLES)):
			shutil.copyfile(SAMPLES[i]+"/"+SAMPLES[i]+'_screen.png','Screen/'+SAMPLES[i]+'_screen.png')				
		shell("cd QcImages/; Rscript ../scripts/RFastQPlot.r && echo 'RFastQPlot.r worked' || echo 'RFastQPlot.r failed'; cd ../Screen/; Rscript ../scripts/RScreenPlot.r  && echo 'RScreenPlotnPlot.r worked' || echo 'RScreenPlot.r failed'; cd ../") 