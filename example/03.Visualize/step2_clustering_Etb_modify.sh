perl /panfs/home/kang/yangpc//bin/EnrichPipeline/EACO_r20150101/bin/run_WGCNA.pl --expr /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/enrich.gsea.pv.tb.logp --expCut 0 --TOMType signed  --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna --step 12  --beta 1

perl /panfs/home/kang/yangpc//bin/EnrichPipeline/EACO_r20150101/bin/cluster_based_on_numeric_matrix.pl --mat /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/enrich.gsea.pv.tb.logp --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster/ --TOMdis /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.networkConstruction.TOM.RData --test


