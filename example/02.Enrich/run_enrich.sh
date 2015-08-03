perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/extract_gene_level_statistics_from_gene_diff.pl --gdiff /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/lmi.gene_exp.diff --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/

perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/02.eaco.enrich.pl  --gmt /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt  --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/ --mcls GSEA --method gsea --rnklist /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/gstat.list --minG 4 --maxG 3000

perl ~/bin/EnrichPipeline/EACO_r20150201/bin/extract_DEG_ids_from_gene_diff.pl --diff /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/lmi.gene_exp.diff  --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/  --padj 0.01

ls /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/*[np] |perl -ne 'chomp;($id)=$_=~/stat\/(\S+)/;print "$id\t$_\n"' > 00.gstat/degid.list

perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/02.eaco.enrich.pl --gmt /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt  --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/ --mcls SGA --method FisherChiSquare --minG 4 --maxG 3000 --idlist /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/degid.list 																																											 
ls /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/SGA/*difsga > SGA/difsga.list

perl /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/bin/cat_enrich_result_to_pvtb.pl --difcfg /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/SGA/difsga.list --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/SGA/01.merge

