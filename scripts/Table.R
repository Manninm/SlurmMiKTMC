#!/usr/bin/env Rscript

library('jsonlite')
mytable<-fromJSON('samples.json')
data<-as.data.frame(mytable)
data <- lapply(mytable, function(x) {
x[sapply(x, is.null)] <- NA
unlist(x)
})
data<-do.call("rbind",data)
Samples<-as.vector(row.names(data))
Directory<-as.vector(paste(Samples,'/',sep=''))
Disease<-rep('Sample',length(Samples))
Batch<-rep('Batch',length(Samples))
groups<-cbind(Directory,Samples,Disease,Batch)
write.table(groups,'groups.txt',quote=FALSE, sep='\t',row.names=FALSE)