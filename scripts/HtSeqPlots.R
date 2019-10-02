#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
groups<-read.table('groups.txt',header=TRUE,sep='\t')
exp<-read.table('HtSeqCounts/CountsGt0_voom.txt',header=TRUE,sep='\t',row.names=1,check.names=FALSE)
names(exp)<-gsub("X","",names(exp))
names(exp)<-gsub('FPKM\\.',"",names(exp))
names(exp)<-gsub("htseq_nopos.inter.str.txt","",names(exp))
names(exp)<-gsub("\\_htseqout_noposuniontxt","",names(exp))
names(exp)<-gsub('FPKM.Sample_','',names(exp))
names(exp)<-gsub('Sample_','',names(exp))
names(exp)<-gsub('*_merged',"",names(exp))
names(exp)<-gsub("\\.","-",names(exp))
names(exp)<-gsub("X","",names(exp))
#make boxplot
length(which(rowSums(exp)>0))
boxexp<-exp[which(rowSums(exp)>0),]
png(paste("Plots/BoxPlotHtSeqVoomGt0.png",width=2000,height=1000)
par(cex.axis=0.9,mai=c(2.8,0.82,0.82,0.42))
boxplot(boxexp,las=2,main=paste("BoxploxNormalizedLog2expression"))
dev.off()
#make pca
names(groups)<-c('Sample','Disease',"Batch")
groups$Sample<-gsub('*_merged','',groups$Sample)
groups$Sample<-gsub('Sample_','',groups$Sample)
exp[is.na(exp)] <- 0 
tmp<-match(colnames(exp),groups[,1])
if(any(is.na(tmp))){return('Colnames of Expression matrix do not match groupfile samples')}
exp<-exp[,(match(groups$Sample,names(exp)))]
exp.pca<-prcomp(t(exp))
labels<-c(as.vector(groups$Disease))
labels<-gsub('ref','REF',labels)
labels<-gsub('Ref','REF',labels)
na<-c(as.character(groups$Sample))
size<-c(as.character(groups$Batch))
Dis<-(unique(as.vector(labels[!(labels=='REF')])))
LabCol<-rainbow(length(Dis))
labels<-gsub('REF','#000000',labels)
cat('Picking Colors')
print(paste(LabCol))
for (i in 1:length(Dis)){
	print(paste(Dis[i]))
	labels<-gsub(Dis[i],LabCol[i],labels)
}
Dis<-append(Dis,'REF')
LabCol<-append(LabCol,'#000000')
batch<-(unique(size))
shape<-c(sample(0:25,size=length(batch)))
cat('Picking Shapes')
print(paste(shape))
print(length(shape))
print(size)
print(batch)
for (i in 1:length(shape)){
	print(paste(i))
	print(paste(shape[i]))
	print(batch[i])
	size<-gsub(paste(batch[i],'\\>',sep=''),shape[i],size)
}
print(paste(size))
pdf(paste("PCA1v2&2v3&3v4_HtSeqVoom_filtered.pdf", width=15, height=15)
plot(exp.pca$x[,1], exp.pca$x[,2], col=labels,pch=c(as.numeric(size)), main="PCA1vs2RefHtSeqVoomFiltered", xlab = "PCA 1", ylab = "PCA 2")
text(exp.pca$x[,1], exp.pca$x[,2], labels=na, pos= 3) #labels points
legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
plot(exp.pca$x[,2], exp.pca$x[,3], col=labels, pch=c(as.numeric(size)), main="PCA1vs2RefHtSeqVoomFiltered", xlab = "PCA 2", ylab = "PCA 3")
text(exp.pca$x[,2], exp.pca$x[,3], labels=na, pos= 3) #labels points
legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
plot(exp.pca$x[,3], exp.pca$x[,4], col=labels, pch=c(as.numeric(size)), main="PCA1vs2RefHtSeqVoomFiltered", xlab = "PCA 3", ylab = "PCA 4")
text(exp.pca$x[,3], exp.pca$x[,4], labels=na, pos= 3) #labels points
legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
dev.off()
#make cluster
noProcs<-args[1]
clustering<-"ward"	
dist<-"euclidian"
bootstrapping=args[2]

library(snow)
library(pvclust)
library(tools)	
library(corrplot)
dim(exp)
length(which(rowSums(exp)>0))
read.table('HtSeqCounts/CountsGt0_voom.txt',header=TRUE,sep='\t',row.names=1,check.names=FALSE)

cl<-makeCluster(noProcs,type="SOCK")
res<-parPvclust(cl,exp,method.hclust="ward.D", method.dist=dist,nboot=bootstrapping)
pdf(paste("HtSeqVoomFilted_ward.D_",bootstrapping,".pdf",sep=""), width=30, height=20)
plot(res)
dev.off()

res<-parPvclust(cl,exp,method.hclust="ward.D2", method.dist=dist,nboot=bootstrapping)
pdf(paste("HtSeqVoomFilted_Ward.D2_",bootstrapping,".pdf",sep=""), width=30, height=20)
plot(res)
dev.off()

stopCluster(cl)

co<-cor(exp)
pdf(paste("Corrplot_HtSeqVoomFilted.pdf",sep=""), width=40, height=40)
corrplot.mixed(co,upper="circle", lower="number", order="hclust",tl.pos="lt")
dev.off()