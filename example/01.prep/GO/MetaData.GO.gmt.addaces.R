

library(GO.db)

bp.aces <- as.list(GOBPANCESTOR)
cc.aces <- as.list(GOCCANCESTOR)
mf.aces <- as.list(GOMFANCESTOR)
all.aces <- c(bp.aces,cc.aces,mf.aces)

bp.off <- as.list(GOBPOFFSPRING)
cc.off <- as.list(GOCCOFFSPRING)
mf.off <- as.list(GOMFOFFSPRING)
all.off <- c(bp.off,cc.off,mf.off)

goterm <- as.list(GOTERM)

tmp <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/00.data/MetaData.GO.gmt",sep="\n",what="character")
go2g <- sapply(tmp,function(x){
  x1 <- strsplit(x,"\t")[[1]]
  x1[3:length(x1)]
})
names(go2g) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
go2g <- go2g[intersect(names(go2g),names(goterm))]

go2g.aces <- setdiff(unique(as.character(unlist(all.aces[names(go2g)]))),names(go2g))
go2g.aces.list <- sapply(go2g.aces,function(x){
		unique(as.character(unlist(go2g[intersect(all.off[[x]],names(go2g))])))
})
go2g <- c(go2g,go2g.aces.list)
go2g <- go2g[setdiff(names(go2g),c("all","GO:0003674","GO:0008150","GO:0005575"))]

#output the gmt file
gmt <- data.frame(GO=names(go2g),
		Term=as.character(sapply(names(go2g),function(x) Term(goterm[[x]]))),
		genes=as.character(sapply(go2g,function(x) paste(x,collapse="\t"))))
write.table(gmt,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO//MetaData.GO.gmt.addaces.gmt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

res <- c(length(go2g),length(unique(as.character(unlist(go2g)))))
names(res) <- c("totalGenSets","totalGenNum")
write.table(res,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO//MetaData.GO.gmt.addaces.gmt.stat",sep="\t",quote=FALSE)

q('no')
