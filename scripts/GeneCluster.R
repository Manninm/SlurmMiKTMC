#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
file<-args[1]
exp<-read.table(file,header=TRUE,sep='\t',row.names=1,check.names=FALSE)
outname<-gsub('.txt', '',file,)
outfile<-paste('../Plots/',outname,sep='')
noProcs<-args[2]
print(noProcs)
clustering<-"ward"	
dist<-"euclidian"
bootstrapping=args[3]
print(bootstrapping)
print(outfile)

library(snow)
library(pvclust)
library(tools)	
library(corrplot)
dim(exp)
length(which(rowSums(exp)>0))

#cl<-makeCluster(noProcs,type="SOCK")
res<-pvclust(exp,method.hclust='ward.D',parallel=as.integer(noProcs),method.dist=dist,nboot=as.numeric(bootstrapping))
pdf(paste(outfile,".Ward.D.",bootstrapping,".pdf",sep=""), width=30, height=20)
plot(res)
dev.off()

res<-pvclust(exp,method.hclust="ward.D2",parallel=as.integer(noProcs),method.dist=dist,nboot=as.numeric(bootstrapping))
pdf(paste(outfile,".Ward.D2.",bootstrapping,".pdf",sep=""), width=30, height=20)
plot(res)
dev.off()



co<-cor(exp)
pdf(paste(outfile,"CorrPlot.pdf",sep=""), width=40, height=40)
corrplot.mixed(co,upper="circle", lower="number", order="hclust",tl.pos="lt")
dev.off()