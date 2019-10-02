library(EBImage)
graphs<-c("adapter_content","per_base_quality", "per_base_sequence_content", "duplication_levels", "per_sequence_gc_content", "kmer_profiles")


graphs<-c("adapter_content","per_base_quality", "per_base_sequence_content", "duplication_levels", "per_sequence_gc_content", "kmer_profiles")

for (graph in graphs) {
	foo<-list()

	filenames<-dir(pattern=paste("*",graph,"*",".png",sep=""))
	for(j in 1:length(filenames)) foo[[j]]<-readImage(filenames[j],type='png')
	pdf(paste(graph,"_pngs.pdf",sep=""),width=60, height=40)
	#pdf("TEST_pngs.pdf", width=60, height=40)
	par(mfrow=c(4,4),mar=c(1,0,8,0))			# mar() is bottom, left, top, right
	for (j in 1:length(filenames)) {
		display(foo[[j]],method="raster" )
		text(x=200, y=0,adj=c(0,1), label=filenames[[j]],cex=2)	#adj() does tome x, y adjustment of the label
	}
	dev.off()
	
}
	