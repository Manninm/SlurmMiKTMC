args = commandArgs(trailingOnly=TRUE)
file<-args[1]
exp<-read.table(file,header=TRUE,sep='\t',row.names=2,check.names=FALSE)
exp<-exp[,-c(1,2)] #For Transcript level
outname<-gsub('.txt', '',file,)
outfile<-paste('../Plots/',outname,sep='')
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
png(paste(outfile,"BoxPlot.png",sep=''),width=2000,height=1000)
par(cex.axis=0.9,mai=c(2.8,0.82,0.82,0.42))
boxplot(boxexp,las=2,main=paste(outname,"BoxploxNormalizedLog2expression",sep=''))
dev.off()