
options(stringsAsFactors=FALSE)
tmp <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO/MetaData.GO.gmt.addaces.gmt",sep="\n",what="character")
set2g <- sapply(tmp,function(x){
		x1 <- strsplit(x,"\t")[[1]]
		x1[3:length(x1)]
		})
names(set2g) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])

set2g.len <- sapply(set2g,length)
brks <- seq(5,max(set2g.len)+5,5)
dt <- matrix(0,nr=length(brks),nc=5)
rownames(dt) <- brks
colnames(dt) <- c("brks","SetNum","SetNumPerc","GeneNum","GeneNumPerc")
dt[,"brks"] <- brks
dt[,"SetNum"] <- as.numeric(sapply(brks,function(x){
			length(set2g.len[set2g.len<x])
		}))
dt[,"SetNumPerc"] <- 100*dt[,"SetNum"]/length(set2g.len)
dt[,"GeneNum"] <- as.numeric(sapply(brks,function(x){
			length(unique(as.character(unlist(set2g[set2g.len<x]))))
			}))
dt[,"GeneNumPerc"] <- 100*dt[,"GeneNum"]/length(unique(as.character(unlist(set2g))))

write.table(dt,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO//MetaData.GO.gmt.addaces.gmt.maxG_and_distinctGeneNum.tb",sep="\t",row.names=FALSE)
dt.long <- data.frame(brks=rep(dt[,"brks"],2),
		Percentage=as.numeric(dt[,c("SetNumPerc","GeneNumPerc")]),
		Group=rep(c("SetNumPerc","GeneNumPerc"),each=NROW(dt))
		)

save(dt.long,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO//MetaData.GO.gmt.addaces.gmt.data.RData")
library(ggplot2)
p.line <- ggplot(dt.long,aes(brks,Percentage,colour=Group))+geom_line()+
	scale_x_log10(breaks=c(1,5,10,100,500,800,1000,1500,2000,2500,3000,4000,5000,10000))+
	xlab("GeneNumber cutoff")+xlim(0,max(brks))+
	theme(axis.text.x = element_text(angle=30))

pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep//GO//MetaData.GO.gmt.addaces.gmt.maxG_and_distinctGeneNum.tb.plot.pdf")
plot(p.line)
dev.off()

q('no')
