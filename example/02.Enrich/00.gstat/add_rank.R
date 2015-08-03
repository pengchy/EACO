
glist = read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/02.Enrich/00.gstat//gstat.list",sep="\t",stringsAsFactors=FALSE)
for(i in glist[,2]){
	a = read.table(i,stringsAsFactors=FALSE)
	a.res = data.frame(a,rank(a[,4],ties.method = "random"))
	i.rnk = paste(i,".rnk",sep="")
	write.table(a.res,file=i.rnk,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	file.copy(from=i.rnk,to=i,overwrite = TRUE)
	file.remove(i.rnk)
	write.table(a.res[,c(1,5)],file=i.rnk,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
