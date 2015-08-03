
date()
source("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/EnrichSGA.R")
source("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/sort.data.frame.R")
source("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/Fisher.Chi.test.R")
supplyID <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/Gegg-VS-Segg.down",what="character",sep="\n")
if(identical("NullFile","NullFile")){
	res.dt <- EnrichSGA(gmt="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt",supplyID=supplyID,p.adjust.methods="fdr",
			test.method="FisherChiSquare",enrichFile="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//SGA///Gegg-VS-Segg.down.difsga")
}else{
	univerID <- scan("NullFile",what="character",sep="\n")
	res.dt <- EnrichSGA(gmt="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt",supplyID=supplyID,univerID=univerID,p.adjust.methods="fdr",
			test.method="FisherChiSquare",enrichFile="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//SGA///Gegg-VS-Segg.down.difsga")
}
date()
q('no')

