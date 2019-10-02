#This Script will output a PCA for a given gene expression matrix with genes as rows, and samples as columns with an accompanied pheno table, which can consist of <= 4 columns consisting of a directory column (which is used for ballgown), sample-name column, category/disease column, and batch/cohort column. Samples in the Sample Column must match colunmn names in Gene expression matrix. The Script will attempt to clean for strings commonly seen in either HtSeq or Stringtie. Given Batch=FALSE, a pca will be colored by Disease/Condition. If Batch=TRUE, then shapes will be used to identify Batch. The Tissue/Feature options are for easily nameing multiple files via conditions and tissue. 
PcaPlot<-function(ExpMat,GroupFile,transcript=FALSE,batch=FALSE,log=FALSE,Directory,features='ExperimentalFeatures',tissue='TissueUsedInAnalysis'){
	#ExpMat='TabDelimited Gene Matrix with Genes as rows/Samples as Columns';GroupFile='Phenotype table <= four columns titles Directory, Sample, Disease, Batch';Batch='If true, will use batch column to assign shapes to different batches, must be < 25 batchs';transcript='use if matrix has both geneids and transcript ids assumes txIDs are 2nd column';Directory='use true if your first column consists of directory paths to samples, often used in Ballgown';features='Naming convenience';tissue='Naming Convenience';Log='use if you wish to log2(+1) transform your data.' 
	library(ggfortify)
	if(Directory){
		groups<-read.table(GroupFile,header=TRUE)
		groups<-groups[,-c(1)]
		groups$Disease<-as.vector(groups$Disease)
	}
	else{
		groups<-read.table(GroupFile,header=TRUE)
	}  
	if(transcript){
		exp<-read.table(ExpMat,header=TRUE,row.names=2)
		exp<-exp[,-c(1,2)]
		names(exp)<-gsub("X","",names(exp))
		names(exp)<-gsub('FPKM\\.',"",names(exp))
		names(exp)<-gsub("htseq_nopos.inter.str.txt","",names(exp))
		names(exp)<-gsub("\\_htseqout_noposuniontxt","",names(exp))
		names(exp)<-gsub('FPKM.Sample_','',names(exp))
		names(exp)<-gsub('Sample_','',names(exp))
		names(exp)<-gsub('*_merged',"",names(exp))
		names(exp)<-gsub("\\.","-",names(exp))
		names(exp)<-gsub("X","",names(exp))
	}
	else{
		exp<-read.table(ExpMat,header=TRUE,row.names=1)
		names(exp)<-gsub("X","",names(exp))
		names(exp)<-gsub('FPKM\\.',"",names(exp))
		names(exp)<-gsub("htseq_nopos.inter.str.txt","",names(exp))
		names(exp)<-gsub("\\_htseqout_noposuniontxt","",names(exp))
		names(exp)<-gsub('FPKM.Sample_','',names(exp))
		names(exp)<-gsub('Sample_','',names(exp))
		names(exp)<-gsub('*_merged',"",names(exp))
		names(exp)<-gsub("\\.","-",names(exp))
		names(exp)<-gsub("X","",names(exp))
	}
	if(log){
		exp<-log2(exp+1)
		}
	else{
	warning('Proceeding without Log Transformation')
	}
	if(batch){
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
		pdf(paste("PCA1v2&2v3&3v4_",tissue,features,".pdf",sep=""), width=15, height=15)
		plot(exp.pca$x[,1], exp.pca$x[,2], col=labels,pch=c(as.numeric(size)), main=paste("PCA1vs2SRC_", features, tissue,sep=""), xlab = "PCA 1", ylab = "PCA 2")
		text(exp.pca$x[,1], exp.pca$x[,2], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
		plot(exp.pca$x[,2], exp.pca$x[,3], col=labels, pch=c(as.numeric(size)), main=paste("PCA2vs3_SRC", features, tissue,sep=""), xlab = "PCA 2", ylab = "PCA 3")
		text(exp.pca$x[,2], exp.pca$x[,3], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
		plot(exp.pca$x[,3], exp.pca$x[,4], col=labels, pch=c(as.numeric(size)), main=paste("PCA3vs4_SRC", features, tissue,sep=""), xlab = "PCA 3", ylab = "PCA 4")
		text(exp.pca$x[,3], exp.pca$x[,4], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		legend("bottomleft",legend=paste(batch,sep=''),pch=as.numeric(shape))
		dev.off()
	
	}
	else{
		names(groups)<-c('Sample','Disease')
		groups$Sample<-gsub('*_merged','',groups$Sample)
		groups$Sample<-gsub('Sample_','',groups$Sample)
		exp[is.na(exp)] <- 0 
		tmp<-match(colnames(exp),groups[,1])
		if(any(is.na(tmp))){return('Colnames of Expression matrix do not match groupfile samples')}
		exp<-exp[,(match(groups$Sample,names(exp)))]
		exp.pca<-prcomp(t(exp))
		labels<-c(as.character(groups$Disease))
		labels<-gsub('ref','REF',labels)
		labels<-gsub('Ref','REF',labels)
		na<-c(as.character(groups$Sample))
		Dis<-(unique(as.vector(labels[!(labels=='REF')])))
		labels<-gsub('REF','#000000',labels)
		LabCol<-rainbow(length(Dis))
		cat('Picking Colors')
		for (i in 1:length(Dis)){
			print(paste(i))
			labels<-gsub(Dis[i],LabCol[i],labels)
		}
		Dis<-append(Dis,'REF')
		LabCol<-append(LabCol,'#000000')
		pdf(paste("PCA1v2&2v3v3&4_",tissue,"_",features,".pdf",sep=""), width=20, height=20)
		plot(exp.pca$x[,1], exp.pca$x[,2], col=labels, pch=16, main=paste("PCA1vs2_", features, tissue,sep=""), xlab = "PCA 1", ylab = "PCA 2")
		text(exp.pca$x[,1], exp.pca$x[,2], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		plot(exp.pca$x[,2], exp.pca$x[,3], col=labels, pch=16, main=paste("PCA2vs3_", features, tissue,sep=""), xlab = "PCA 2", ylab = "PCA 3")
		text(exp.pca$x[,2], exp.pca$x[,3], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		plot(exp.pca$x[,3], exp.pca$x[,4], col=labels, pch=16, main=paste("PCA3vs4_", features, tissue,sep=""), xlab = "PCA 3", ylab = "PCA 4")
		text(exp.pca$x[,3], exp.pca$x[,4], labels=na, pos= 3) #labels points
		legend("bottomright",legend=paste(Dis,sep=''),fill=paste(LabCol,sep=''))
		dev.off()
		
	}
}  