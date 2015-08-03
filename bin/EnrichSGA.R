EnrichSGA <- function(gmt,supplyID,univerID,
                       p.adjust.methods=c("holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                       enrichFile="enrichOutfile.difsga",test.method=c("FisherChiSquare","HyperGeometric"),
                       filt=c("p","adjp"),pc=0.05,Unanno=c("Y","N"))
{
#set name to description
	tmp <- scan(gmt,sep="\n",what="character")
	gSets <- sapply(tmp,function(x){
			x1 <- strsplit(x,"\t")[[1]]
			x1[3:length(x1)]
		})
	names(gSets) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
	setDes <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][2])
	names(setDes) <- names(gSets)

  ##set to genes and gene to sets
  gene2sets <- tapply(rep(names(gSets),sapply(gSets,length)),
                     as.character(unlist(gSets)),unique)
  
  if(missing(p.adjust.methods)) p.adjust.methods="fdr"
  if(missing(univerID)) univerID <- names(gene2sets)
  if(missing(filt)) filt <- "adjp"
	if(missing(Unanno) | identical(Unanno,"N")){
  	univerID <- intersect(univerID,names(gene2sets))
 	 ids <- intersect(supplyID,univerID)
 	 totAnnNum <- length(univerID)
  	annoDiffNum <- length(ids)
	}else{
		totAnnNum <- length(univerID)
		annoDiffNum <- length(supplyID)
  	univerID <- intersect(univerID,names(gene2sets))
 	 ids <- intersect(supplyID,univerID)
	}
  if(length(univerID) != length(gene2sets)){
    map2id <- lapply(gSets,function(x) intersect(x,univerID))
    map2id <- map2id[sapply(map2id,length)>0]
  }else{
    map2id <- gSets
  }
  mapNum <- length(map2id)
  diffSet.tbl <-
    data.frame("gSetID"=names(map2id),
               "gSetDes"=as.character(setDes[names(map2id)]),
               "Pvalue"=rep(0,mapNum),
							 "AdjustedPv"=rep(0,mapNum),
               "x"=rep(0,mapNum),
               "y"=rep(0,mapNum),
               "n"=rep(0,mapNum),
               "N"=rep(0,mapNum),
               "EnrichDirect"=rep("a",mapNum),
               "GeneIDs"=rep("a",mapNum),stringsAsFactors=FALSE)
  rownames(diffSet.tbl) <- names(map2id)
  for(j in rownames(diffSet.tbl)){
    knum <- length(intersect(ids,map2id[[j]]))
    if(knum <1) next
    dif.data <- c(knum,length(map2id[[j]]),annoDiffNum,totAnnNum)
    diffSet.tbl[j,5:8] <- dif.data
    knum.rm <- length(map2id[[j]])-knum
    tot.rm <- totAnnNum-annoDiffNum-knum.rm
    dif.data <- c(knum,knum.rm,annoDiffNum-knum,tot.rm)
    mt <- matrix(dif.data,nr=2)
    if(missing(test.method)) test.method="FisherChiSquare"
    test.pv <- switch(test.method,
                      FisherChiSquare=Fisher.Chi.test(mt)$p.v,
                      HyperGeometric=hyper.test(mt)
                      )
    diffSet.tbl[j,3] <- test.pv
    diffSet.tbl[j,9] <- ifelse((dif.data[1]/dif.data[2])>(dif.data[3]/dif.data[4]),
                                "Over","Under")
    diffSet.tbl[j,10] <- paste(intersect(ids,map2id[[j]]),collapse=",")
  }
  diffSet.tbl <- diffSet.tbl[diffSet.tbl[,5]>0,]
  diffSet.tbl[,4] <- p.adjust(diffSet.tbl[,3],method=p.adjust.methods)
  diffSet.tbl <- sort.data.frame(diffSet.tbl,key="Pvalue")
  write.table(diffSet.tbl,file=enrichFile,sep="\t",quote=FALSE,row.names=FALSE)
  enrichFile.filt <- paste(enrichFile,".filt",sep="")
  if(filt=="adjp")
    diffSet.filt <- diffSet.tbl[diffSet.tbl[,"AdjustedPv"] < pc & diffSet.tbl[,"EnrichDirect"]=="Over",]
  else
    diffSet.filt <- diffSet.tbl[diffSet.tbl[,"Pvalue"] < pc & diffSet.tbl[,"EnrichDirect"]=="Over",]  
  write.table(diffSet.filt,file=enrichFile.filt,sep="\t",quote=FALSE,row.names=FALSE)
  return(diffSet.tbl)
}
