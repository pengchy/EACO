#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

Zscore <- function(n,r,N,R){
	(r-n*R/N)/sqrt(n*(R/N)*(1-R/N)*(1-(n-1)/(N-1)))
}

ZGcomb <- function(n,r,N,R){
	Zscore(n,r,N,R) * log2(r)
}

read.gmt <- function(file){
	tmp <- scan(file,what="character",sep="\n")
	gsets <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][-c(1,2)])
	names(gsets) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
	setDes <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][2])
	names(setDes) <- toupper(names(gsets))
	names(gsets) <- toupper(names(gsets))
	list(gsets=gsets,setDes=setDes)
}

write.gmt <- function(gsets,setDes,file){
	gmt.dt <- data.frame(names(gsets),
			setDes,
			sapply(gsets,function(x) paste(x,collapse="\t"))
			)
	write.table(gmt.dt,file=file,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
}


GOFunction.local <- function(
		local.tb, #in the format: ancesID offID ancesGnum offGnum spcNum
		orig.enrich, #the original enrichment result table, used to filter
		sub.enrich, #the sub gene sets enrichment result table, every column name should be same to the orig.enrich
		orig.gmt, #gmt file corresponding to orig.enrich
		pv=0.05, #p value cut off to filter the orig.enrich and sub.enrich
		out.dir="./" #output directory
		){
	delid <- vector("character")
	goids <- unique(as.character(c(rownames(local.tb),unlist(sapply(local.tb[,1],function(x) strsplit(x,",")[[1]])))))
	replid <- vector("list",length=length(goids))
	names(replid) <- goids

	for(i in rownames(local.tb)){
		offids <- strsplit(local.tb[i,1],",")[[1]]
		mergids <- c(i,offids)

#delete ancestor if ancestor have less than 3 specific genes
		if(local.tb[i,4]<3 || 
				(local.tb[i,"spcNum"]>=3 && !is.element(i,rownames(sub.enrich)))){
			delid <- c(delid,i)
			for(j in offids){
				replid[[j]] <- append(replid[[j]],i)
			}
			next
		}
		
		if(local.tb[i,"spcNum"]>=3 && !is.element(i,rownames(sub.enrich))){
			cat("this term may be filtered by piano because the minimal gene numer limitation or these specfici genes not included in the analyzed genes",i,"\n")
		}

#delete ancestor
		if(length(which(apply(orig.enrich[mergids,],2,function(x) all(x<pv)))) == #at least one sample satisfy the pvalue criteria
	  	 length(which(apply(orig.enrich[mergids,],2,function(x) any(x<pv))))){
			pos <- which(apply(orig.enrich[mergids,],2,function(x) all(x<pv)))
			pos.judge <- sapply(pos,function(x){
				all(sub.enrich[i,pos]>pv)
				})
			
			if(all(pos.judge)){
				delid <- c(delid,i)
				for(j in offids){
					replid[[j]] <- append(replid[[j]],i)
				}
			}
		}

#delete offspring
		pos <- which(apply(orig.enrich[mergids,],2,function(x) any(x<pv)))
		pos.judge <- sapply(pos,function(x) {
				all(orig.enrich[i,x]<orig.enrich[offids,x])
			})
		if(all(pos.judge) && all(sub.enrich[i,pos]<pv)){
			delid <- c(delid,offids)
			replid[[i]] <- append(replid[[i]],offids)
		}
	}

	delid <- unique(delid)
	replid <- sapply(replid, unique)
	replid <- replid[sapply(replid,length)>0]
	replid.dt <- data.frame(
			    remain=rep(names(replid),sapply(replid,length)),
					    del=as.character(unlist(replid)))
	replid <- replid[setdiff(names(replid),delid)]

	#output filtered gmt file
	gsets.afterlocal <- orig.gmt[["gsets"]][setdiff(names(orig.gmt[["gsets"]]),delid)]
	gsets.afterlocal <- gsets.afterlocal[intersect(rownames(orig.enrich),names(gsets.afterlocal))]
	setDes.afterlocal <- orig.gmt[["setDes"]][setdiff(names(orig.gmt[["setDes"]]),delid)]
	setDes.afterlocal <- setDes.afterlocal[intersect(rownames(orig.enrich),names(setDes.afterlocal))]
	setDes.afterlocal[names(replid)] <- sapply(names(replid),function(x){
			    x1 <- paste(replid[[x]],orig.gmt[["setDes"]][replid[[x]]],sep=":",collapse="|")
			    x2 <- paste(orig.gmt[["setDes"]][x],x1,sep="|")
			    x2})
	write.gmt(gsets.afterlocal,setDes.afterlocal,file=paste(out.dir,"/afterlocal.filt.gmt",sep=""))
	res <- list()
	res[["gsets"]] <- gsets.afterlocal
	res[["setDes"]] <- setDes.afterlocal
	res[["replid"]] <- replid
	return(res)
}


GOFunction.glb <- function(glb.tb,orig.enrich,sub.enrich,orig.gmt,pv=0.05,out.dir="./"){
	delid <- vector("character")
	goids <- unique(c(glb.tb[,1],glb.tb[,2]))
	replid <- vector("list",length=length(goids))
	names(replid) <- goids
	for(i in 1:dim(glb.tb)[1]){
		twoids <- as.character(glb.tb[i,1:2])
		threeN <- as.numeric(glb.tb[i,3:5])
 #only treat the set pairs with their specific gene number >= 3
		if(threeN[2] < 3 || threeN[3] < 3) next
#only treat the set pairs both enriched at one sample
		pos1 <- as.numeric(which(apply(orig.enrich[twoids,],2,function(x) all(x<pv))))
		pos2 <- as.numeric(which(apply(orig.enrich[twoids,],2,function(x) any(x<pv))))
		if(length(pos1) != length(pos2)) next 
		if(length(pos1)==0) cat("find 0 position\n")

#disregard two sets with the difference between two specific gene number greater than 500
		if(abs(threeN[2]-threeN[3])>500) next

		twoids.spc <- twoids
		twoids.spc[1] <- paste(twoids[1],"_SPC_VS_",twoids[2],sep="")
		twoids.spc[2] <- paste(twoids[1],"_VS_",twoids[2],"_SPC",sep="")

#select the sets with more genes as representative when their specific gene were both not be included by the analyzed gene list
		if(!is.element(twoids.spc[1],rownames(sub.enrich)) && !is.element(twoids.spc[2],rownames(sub.enrich))){
			cat("this two sets may be have small specific gene and not included in the analyzed gene list: ",twoids.spc,"\n")
			pos3 <- which.max(sapply(orig.gmt[["gsets"]][twoids],length))
			delid <- append(delid,twoids[pos3])
			replid[[twoids[-pos3]]] <- append(replid[[twoids[-pos3]]],twoids[pos3])
			next
		}
#or one set's specific gene not be included by the analyzed gene list
		if(!is.element(twoids.spc[1],rownames(sub.enrich)) && is.element(twoids.spc[2],rownames(sub.enrich))){
			delid <- append(delid,twoids[1])
			replid[[twoids[2]]] <- append(replid[[twoids[2]]],twoids[1])
			next
		}
		if(!is.element(twoids.spc[2],rownames(sub.enrich)) && is.element(twoids.spc[1],rownames(sub.enrich))){
			delid <- append(delid,twoids[2])
			replid[[twoids[1]]] <- append(replid[[twoids[1]]],twoids[2])
			next
		}

		if(all(apply(data.frame(sub.enrich[twoids.spc,pos1]),2,function(x) x[1] > pv && x[2] < pv))){
			delid <- append(delid,twoids[1])
			replid[[twoids[2]]] <- append(replid[[twoids[2]]],twoids[1])
		}
	
		if(all(apply(data.frame(sub.enrich[twoids.spc,pos1]),2,function(x) x[2] > pv && x[1] < pv))){
			delid <- append(delid,twoids[2])
			replid[[twoids[1]]] <- append(replid[[twoids[1]]],twoids[2])
		}
	}

	delid <- unique(delid)
	replid <- sapply(replid, unique)
	replid <- replid[sapply(replid,length)>0]
	replid <- replid[setdiff(names(replid),delid)]
	
#output filtered gmt file
	gsets.afterlocal <- orig.gmt[["gsets"]][setdiff(names(orig.gmt[["gsets"]]),delid)]
	gsets.afterlocal <- gsets.afterlocal[intersect(rownames(orig.enrich),names(gsets.afterlocal))]
	setDes.afterlocal <- orig.gmt[["setDes"]][setdiff(names(orig.gmt[["setDes"]]),delid)]
	setDes.afterlocal <- setDes.afterlocal[intersect(rownames(orig.enrich),names(setDes.afterlocal))]
	setDes.afterlocal[names(replid)] <- sapply(names(replid),function(x){
			x1 <- paste(replid[[x]],orig.gmt[["setDes"]][replid[[x]]],sep=":",collapse="|")
			x2 <- paste(orig.gmt[["setDes"]][x],x1,sep="|")
			x2})
	write.gmt(gsets.afterlocal,setDes.afterlocal,file=paste(out.dir,"/afterGlobal.filt.gmt",sep=""))
	res <- list()
	res[["gsets"]] <- gsets.afterlocal
	res[["setDes"]] <- setDes.afterlocal
	res[["replid"]] <- replid
	return(res)
}




