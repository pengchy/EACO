
date()
#read in the data
tmp <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//Kappa/gSets.gmt.filt.gmt",what="character",sep="\n")
gos <- sapply(tmp,function(x){
		strsplit(x,"\t")[[1]][-c(1,2)]
		})
names(gos) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
gos.comb <- t(combn(names(gos),2))
gos.comb <- rbind(gos.comb,matrix(rep(names(gos),2),nc=2))
geneN <- length(unique(unlist(gos)))
kval <- apply(gos.comb,1,function(x){
		C11 <- length(intersect(gos[[x[1]]],gos[[x[2]]]))
		C00 <- geneN - length(unique(union(gos[[x[1]]],gos[[x[2]]])))
		Cg1 <- length(gos[[x[1]]])
		Cg2 <- length(gos[[x[2]]])
		C10 <- Cg1 - C11
		C01 <- Cg2 - C11
		Oab <- (C11 + C00) / geneN
		Aab <- (Cg1 * Cg2 + (geneN - Cg1) * (geneN - Cg2))/(geneN*geneN)
		Kab <- (Oab - Aab)/ (1- Aab)
		c(Kab,range(c(C11/Cg1,C11/Cg2)))
		})
kval <- round(t(kval),digit=5)

res <- cbind(gos.comb,kval)
colnames(res) <- c("Set1","Set2","Kappa","RatMin","RatMax")
write.table(res,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//Kappa//gSets.gmt.filt.gmt.kappa.tb",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
date()
q('no')

