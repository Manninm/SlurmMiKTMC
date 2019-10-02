# http://r.789695.n4.nabble.com/Loading-an-image-picture-png-jpeg-to-screen-td2244923.html
library(EBImage)

graphs<-c("screen")


for (graph in graphs) {
	foo<-list()
	# comment that out if filled filenames above
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


