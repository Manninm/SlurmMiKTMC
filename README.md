## Greatlakes Issues
This pipeline was established to specifically function on the Greatlakes Server which uses Slurm. There is module system in place on Greatlakes to control user's environments. The module fftw/3.3.8 doesn't play nice with R. fftw/3.3.8, and the library fftwtools is necessary for EBImage to be installed by Bioconductor. There are workarounds for this problem Ken Weiss kindaly made a module REBImage/3.9 that will function with this pipeline. Also, I have included a version of the EBImage source code, with the fftw dependencies removed. **Warning, this is not a fully functional verson of EBImage** Only Image and display are needed for the pipeline, and those are the only functions I tested. You can load this library like anyother if you move it into you local R library directory.
## citation
Much of this pipeline was inspired by https://github.com/snakemake-workflows and https://github.com/crazyhottomy. The fastq2jason.py script was modified from the original by https://github.com/crazyhottomy, but the Snakefile and modularized rules were inspired by https://github.com/snakemake-workflows. All Files in rules and scripts are my own work. If you use this pipeline, please cite Manninm/MiKTMCSnakemakePipeline

## How to use Pipeline
Most of the specifics of the pipeline can be handled in the config.yaml file. The snakefile, rules and cluster.json SHOULD NOT BE EDITED BY HAND. If you absolutley need to edit cluster.json, I recommend https://jsoneditoronline.org/. Snakemake is very sensitive to syntax, and just saving a file in the wrong format can cause problems.

Download the pipeline from Github or transfer the pipeline from my home directory on 76 server

```bash
tar -xvf MiKTMCSnakemakePipeline.tar.gz
mv -v MiKTMCSnakemakePipeline/* .
rm -r MiKTMCSnakemakePipeline/
```

Do dry run to check outputs and rules

```bash
snakmake -npr -s Snakefile
```
Make DAG or Rulegraph
```bash
snakemake --forceall --rulegraph -s Snakefile | dot -Tpng > rulegrap.png
snakemake --forceall --rulegraph -s Snakefile | dot -Tpdf > rulegrap.pdf
snakemake --forceall --dag -s Snakefile | dot -Tpng > dag.png
snakemake --forceall --dag -s Snakefile | dot -Tpdf > dag.pdf
```
Run locally using 22 cores
```bash
snakemake -j 22 -s Snakefile
```
Run on Greatlakes and Slurm
FYI, the --flags used in the snakemake command call must be somewhere in cluster.json, wwether under the default heading, or the rule heading. If --tasks-per-node is called in the command call, and only --tasks-per-cpu is in your default/rule heading, snakemake will complain that "Wildcards have no attribute..."

```bash 
snakemake -j 999 --cluster-config cluster.json --cluster 'sbatch --job-name {cluster.job-name} --ntasks-per-node {cluster.ntasks-per-node} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'
```
You can also what your queued jobs 
```bash
watch squeue -u $(whoami)
```
## The Workflow of Pipeline
The workflow is as seen below
![](rulegraph.png)

The pipeline expects a directory format as the below example
CAUTION Four or more samples must be included, or the PCA scripts will break. It expects pair-end reads. To my knowledge, the pipeline will not accomodate single-end reads.
```bash
RNAseqTutorial/
├── Sample_70160
│   ├── 70160_ATTACTCG-TATAGCCT_S1_L001_R1_001.fastq.gz
│   └── 70160_ATTACTCG-TATAGCCT_S1_L001_R2_001.fastq.gz
├── Sample_70161
│   ├── 70161_TCCGGAGA-ATAGAGGC_S2_L001_R1_001.fastq.gz
│   └── 70161_TCCGGAGA-ATAGAGGC_S2_L001_R2_001.fastq.gz
├── Sample_70162
│   ├── 70162_CGCTCATT-ATAGAGGC_S3_L001_R1_001.fastq.gz
│   └── 70162_CGCTCATT-ATAGAGGC_S3_L001_R2_001.fastq.gz
├── Sample_70166
│   ├── 70166_CTGAAGCT-ATAGAGGC_S7_L001_R1_001.fastq.gz
│   └── 70166_CTGAAGCT-ATAGAGGC_S7_L001_R2_001.fastq.gz
├── scripts
├── groups.txt
└── Snakefile
```
The pipeline uses two types of annotation and feature calling for redundancy in the event that one pipeline fails/gives 'wonky' results
Upon initiating the snakemake file, the snakemake preamble will check fastq file extensions (our lab uses .fq.gz for brevity) and change any fastq.gz to fq.gz. The preamble will then generate a samples.json file using fastq2json.py. You should check samples.json and makesure it is correct because the rest of the pipeline uses this file to create wildcars, which is the driving force behind snakemake.
If no groupfile (groups.txt) was provided, the preample will generate one for you. This file is necessary to run ballgown as well as the PCA plots. This should also be checked for errors. If you provide your own groups.txt, it should be in the format below
```bash
Directory       Samples Disease Batch
Sample_70160/   Sample_70160    Sample  Batch
Sample_70161/   Sample_70161    Sample  Batch
Sample_70162/   Sample_70162    Sample  Batch
Sample_70166/   Sample_70166    Sample  Batch
```
The directory and sample names should correspond and be in the order as they appear in the directory. The sample and batch columns can be used to designate phenotype data and any batchs you may have. If you have varying 'Disease' types, you can then use this file for differential expression and use the batch column to correct for batch affects. The PCA plotting scripts will plot Disease types in different colors, and different Batchs with different shapes

I have attempted to make this pipeline as streamlined and automatic as possible. It could incorporate differential expression, but I feel that the pipeline completes sufficient tasks for review before Differetial Analysis. In the even that a cohort has Glom and Tub samples, it would be wise to run each separately in their own pipeline. Adding another child directory would be more difficult to code rules for. If there are any plots, qc tools or metrics that you use in your personal analysis, those can be integrated upon request.

If you plan on running multiple snakemake pipelines per account/username, you may wish to add a prefix to _default_ job-name in cluster.json, as well as to job-name in Snakemake.sbat, this will allow you to pull specific job-ids from the squeue command. I have an alias in my bash_profile that returns a specific form of squeue
```bash 
alias JobPriority='squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %Q"'squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %Q"
```
In the event that I want to cancel all jobs from a certain pipeline, I can do so by
```bash
all=$(JobPriority | grep 17_ | awk {'print $1'})
scancel $all
```
