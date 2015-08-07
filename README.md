# EACO
Enrichment Analysis for Customized Organism

-----------------
r20150101
-----------------

1. From this version, the EnrichmentPipeline has been adopted by EACO, which is more powerfull and comprehensive.


----------------------------------
Features of EnrichmentPipeline
----------------------------------

The EnrichmentPipeline package has been designed initially to do enrichment
analysis for locust transciptome, which was published on PLoS ONE, 5(12): e15633
at 2010. Because of its easy to use, afterwards, the pipeline was applied to
many projects at BGI.

This pipeline is designed to do enrichment analysis for the category data, 
such as GO/KEGG/IPR, etc. These data are characterised as one class containing
many genes and one gene is involved in many categories. So, for every class,
a p value is calculated representing the probability that the observed numbers
of counts could have resulted from randomly distributing this class between
the tested gene list and the reference gene list. Usually, the reference gene
list is the total genes of one organism annotated, but one can use customized
gene list as reference. The p value can be approximated by ChiSquare-Fisher
test or hypergeometric test. The p values are corrected by fdr or other
methods.

For GO enrichment, first the hierarchy feature of GO information is converted
to MetaGO_release.RData object from GO.db package of bioconductor. This
package will be updated every half year. Second, the GOdata_organism.RData are
constructed for specific organisms. This file contain the map relationship
for GOs at level 3 to genes and genes to the lowest level GO annotated,
several functions used in the enrichment analysis. Finally the enrichment
analysis is done based on the GOdata_organism.RData and supplied gene list at
level 3 and all levels. To remove the redundance, if the GOs enriched at
different levels with parent-child relationship and have the same gene list,
the lowest level is choose and other levels are filtered. The results file is
suffix "merg".

Reference:
*  Chen S, Yang P, Jiang F, Wei Y, Ma Z, Kang L: De Novo Analysis of Transcriptome Dynamics in the Migratory Locust during the Development of Phase Traits. PLoS One 2010, 5(12):e15633.
*  Huang da, W., B.T. Sherman, and R.A. Lempicki. 2009. Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. Nucleic Acids Res 37: 1-13.
*  Beissbarth, T. and T.P. Speed. 2004. GOstat: find statistically overrepresented Gene Ontologies within a group of genes. Bioinformatics 20: 1464-1465.
