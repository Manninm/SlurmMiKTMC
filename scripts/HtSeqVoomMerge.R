#!/usr/bin/env Rscript
file_list<-dir(pattern="_htseq.cnt$", recursive=TRUE)
dataset <- do.call("cbind",lapply(file_list,FUN=function(files){read.delim(files, sep="\t", header=FALSE,row.names=1)}))
file_list<-gsub("_htseq.out_no.pos.union.txt", "", file_list)
colnames(dataset)<-file_list
dim(dataset)
write.table(dataset,"allCounts.txt", sep="\t", quote=FALSE)
length(which(rowSums(dataset)>0))
dataset_short<-dataset[which(rowSums(dataset)>0),]
dataset_short<-dataset_short[-c(grep("__.*",rownames(dataset_short))),]
dim(dataset_short)
write.table(dataset_short, "CountsGt0.txt", sep="\t", quote=FALSE)
library(limma)
exp<-voom(dataset_short)
exp<-2^exp$E
exp<-log2(exp+1)
write.table(exp,"CountsGt0_voom.txt", sep="\t", quote=FALSE)
median<-apply(exp,1,median)
exp<-cbind(exp, median)
exp<-exp[order(median, decreasing=T),]
exp<-subset(exp, select=-c(median))
half<-floor(dim(exp)[1]/2)
low<-exp[half:dim(exp)[1],]
stdev<-sd(as.matrix(low))
minval<-min(low)
cutoff<-minval+(2*stdev)
keep<-apply(exp,1,FUN=function(row){which(length(which(row>=cutoff))>=3)})
filtered<-exp[names(which(keep>0)),]
write.table(filtered, "CountsGt0_voom_filtered.txt", sep="\t", quote=FALSE)
save.image("filteringR.image",compress=T)
