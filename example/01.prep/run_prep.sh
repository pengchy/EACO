perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/01.eaco.prep.pl --GO --gSets /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/MetaData.GO.gmt --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/


cat ../00.data/MetaData.IPR.gmt ../00.data/MetaData.kegg.gmt > gSets.gmt

cat GO/MetaData.GO.gmt.addaces.gmt.filt.gmt >> gSets.gmt

perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/01.eaco.prep.pl --gSets /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/gSets.gmt --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/


perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/create_customized_annotation_package_for_bioconductor.pl --go /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/MetaData.GO.gmt.3cl --genDes /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/gene.id.desc --genChr /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/geneid2scfid --taxid 7004 --genus Locusta --species migratoria --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Org.db


