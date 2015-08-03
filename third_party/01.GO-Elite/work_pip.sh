python ~/soft/05.SystemBiology/GO-Elite_v.1.2.5-Py/GO_Elite.py --species Lm --input  /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/input/ --denom /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/denomitor/ --output /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/output/ --mod EACO --dataToAnalyze KEGG

perl ~/bin/EnrichPipeline/EACO_r20150201/bin/prepare_input_for_GOElite.pl --idlist ~/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/degid.list --outDir /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/input/ --SysAbr En

python ~/soft/05.SystemBiology/GO-Elite_v.1.2.5-Py/GO_Elite.py --species Lm --input  /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/input/ --denom /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/denomitor/ --output /panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/third_party/01.GO-Elite/output/  --dataToAnalyze all --version EnsMart00

