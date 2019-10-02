library(ballgown)
library(RSkittleBrewer)
library(dplyr)
library(devtools)
library(gtools)
dir.create('BallGown')
pheno_data<-read.table("groups.txt",header=TRUE) #table of folder names of samples to be evaluated
#Start Creating Table to subset for Bg Object
dirs<-as.vector(pheno_data[,c(1)])
pheno_data<-pheno_data[,c(2,3,4)]
#Create Table for Disease designations to append to alldata
bg<-ballgown(samples=dirs, pData=pheno_data, meas="all")
gene_expression = gexpr(bg)
write.table(gene_expression, "BallGown/gene_expression_table.txt", quote=FALSE, sep="\t")
transcript_fkm<-texpr(bg,'FPKM')
transcript_fkm<-log2(transcript_fkm+1)
transcript_fkm<-data.frame(geneNames=ballgown::transcriptNames(bg),geneIDs=ballgown::geneIDs(bg),transcript_fkm)
write.table(transcript_fkm, "BallGown/transcript_fpkm.txt",  quote=FALSE, sep="\t")
transcript_cov = texpr(bg, 'cov')
write.table(transcript_cov, "BallGown/transcript_cov.txt", quote=FALSE, sep="\t")
whole_tx_table = texpr(bg, 'all')
write.table(whole_tx_table, "BallGown/whole_tx_table.txt", quote=FALSE, sep="\t")
#new filtering technique
#logtransform2x+1
#extract ballgown object for filtering
exp<-gexpr(bg)
exp<-log2(exp+1)
median<-apply(exp,1,median)
exp<-cbind(exp, median)
# order by median, decreasing
exp<-exp[order(median, decreasing=T),]
exp<-subset(exp, select=-c(median))
#take the lower half of the expression matrix, those have lower expression
half<-floor(dim(exp)[1]/2)
low<-exp[half:dim(exp)[1],]
stdev<-sd(as.matrix(low))
minval<-min(low)
# cutoff is the minimum value +2* stdev of he lower half
# cutoff is the minimum value +2* stdev of he lower half
cutoff<-minval+(2*stdev)
bg_filt<-exprfilter(bg, cutoff=cutoff)
filt_gene_expression = gexpr(bg_filt)
filt_gene_expression<-log2(filt_gene_expression+1)
write.table(filt_gene_expression, "BallGown/filt_gene_expression_table.txt", quote=FALSE, sep="\t")
transcript_fkm<-texpr(bg_filt,'FPKM')
filt_transcript_fkm<-log2(transcript_fkm+1)
filt_transcript_fkm<-data.frame(geneNames=ballgown::transcriptNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),filt_transcript_fkm)
write.table(filt_transcript_fkm, "BallGown/filt_transcript_fpkm.txt", quote=FALSE, sep="\t")
filt_transcript_cov = texpr(bg_filt, 'cov')
write.table(filt_transcript_cov, "BallGown/filt_transcript_cov.txt", quote=FALSE, sep="\t")
filt_whole_tx_table = texpr(bg_filt, 'all')
write.table(filt_whole_tx_table, "BallGown/filt_whole_tx_table.txt", quote=FALSE, sep="\t")
save.image("BallGown/BallGown.Rimage", compress=T)




