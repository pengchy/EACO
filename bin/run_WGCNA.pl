#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--expr    \t<str>\texpression table file with the format, the first line is the header
\t\twithout the first column. every line is one gene
\t\tthis data will be log2 transformed and normalized to variance 1 and mean 0
--exprChek\t<str>\tchecking the sample quality
--expCut  \t<int>\tthe maximum value that the genes expression value smaller than this
\t\tvalue will be filtered (default is 5)
--tp      \t<str>\tthe max power used for test the power threshold (default 20, minimal 12)
--beta    \t<str>\tthe soft thresholding power set according to the Analysis of network
\t\ttopology for various soft-thresholding powers.
--TOMType \t<str>\ta character string specifying TOM type to be calculated. One of 'unsigned',
\t\t\t'signed' (default unsigned)
--merge   \t<str>\twhether or not (Y or N) to merge the modules (default Y)
--outDir  \t<str>\toutput directory
--step    \t<str>\tsteps:
\t\t1. power selection: preprocess expression data and plot the power law graph. then according to the
\t\t\tplot get the power parameter beta and set to the parameter --beta
\t\t2. network construction: using the --beta parameter run network construction and plot clustering dendrogram

\n";

my($expr,$outDir,$step,$beta,$tp,$TOMType,$merge,$expCut,
		$exprChek);

GetOptions(
		"expr:s"=>\$expr,
		"outDir:s"=>\$outDir,
		"step:s"=>\$step,
		"beta:f"=>\$beta,
		"tp:i"=>\$tp,
		"TOMType:s"=>\$TOMType,
		"merge:s"=>\$merge,
		"expCut:i"=>\$expCut,
		"exprChek"=>\$exprChek
);

$step ||= "1";
$tp ||= 20;
$TOMType ||= "unsigned";
die $usage if !defined $expr;
die "not defined beta !!!\n $usage\n" if ($step eq "2" && ! defined $beta);

my $expr_base=basename($expr);
my $rscpt;

`mkdir -p $outDir` if ! -e $outDir;

if($step =~ /1/){
	if(!defined $expCut){
		$expCut ||= 5;
	}
	if(defined $exprChek){
		$exprChek="y";
	}else{
		$exprChek="n";
	}
	$rscpt=qq#
\#reading data
library(WGCNA) 
options(stringsAsFactors = FALSE)
expr <- read.table("$expr",check.names=FALSE)
expr <- expr[apply(expr,1,function(x) any(x>$expCut)),]
\#output the filtered expression value
write.table(expr,file="$outDir/$expr_base.after_filt",sep="\\t",quote=FALSE)
datExpr0 <- as.data.frame(t(expr))
if(identical("y","$exprChek")){
	gsg <- goodSamplesGenes(expr,verbos=3)
	if (!gsg\$allOK){
		if (sum(!gsg\$goodGenes)>0)
			printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg\$goodGenes], collapse = ", ")));
		if (sum(!gsg\$goodSamples)>0)
			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg\$goodSamples], collapse = ", ")));
		datExpr0 = datExpr0[gsg\$goodSamples, gsg\$goodGenes]
	}
}

datExpr0.dimn <- dimnames(datExpr0)
datExpr0 <- log2(datExpr0+1)
datExpr0 <- t(scale(t(datExpr0),center=TRUE,scale=TRUE))
datExpr0 <- scale(datExpr0,center=TRUE,scale=TRUE)
dimnames(datExpr0) <- datExpr0.dimn
\#checking data quality
if(identical("y","$exprChek")){
	sampleTree = flashClust(dist(datExpr0), method = "average");
	pdf("$outDir/$expr_base.Sample_clustering_to_detect_outliers.pdf",height=9,width=12)
	par(cex = 0.6);
	par(mar = c(0,4,2,0))
	plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
		cex.axis = 1.5, cex.main = 2)
	dev.off()
}
datExpr <- datExpr0
\#construct coexpression network
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=$tp, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
\#select beta
pdf("$outDir/$expr_base.Power_selection.pdf",height=5,width=9)
par(mfrow = c(1,2));
cex1 = 0.9; 
\# Scale-free topology fit index as a function of the soft-thresholding power 
plot(sft\$fitIndices[,1], -sign(sft\$fitIndices[,3])*sft\$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft\$fitIndices[,1], -sign(sft\$fitIndices[,3])*sft\$fitIndices[,2],
labels=powers,cex=cex1,col="red");
\# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
\# Mean connectivity as a function of the soft-thresholding power
plot(sft\$fitIndices[,1], sft\$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft\$fitIndices[,1], sft\$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
save(datExpr,expr,sft,file="$outDir/$expr_base.expr.RData")
q('no')
#;
	open O,">","$outDir/$expr_base.expr.R";
	print O $rscpt;
	close O;

	open O,">","$outDir/step1_expr.sh";
	print O "R CMD BATCH $outDir/$expr_base.expr.R $outDir/$expr_base.expr.Rout\n";
	close O;

	`R CMD BATCH $outDir/$expr_base.expr.R $outDir/$expr_base.expr.Rout`;
}

#step2:
#1. construct network based on selected beta
if($step =~ /2/){
	die "undefined beta\n\n" if ! defined $beta;
	$merge ||= "N";
	$rscpt=qq#
library(WGCNA)
load("$outDir/$expr_base.expr.RData")
\#construct network step-by-step
enableWGCNAThreads()
\#Co-expression similarity and adjacency
adjacency <- adjacency(datExpr, power = $beta,type="$TOMType");
\# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dimnames(TOM) <- dimnames(adjacency)
dissTOM = 1-TOM
save(adjacency,dissTOM,file="$outDir/$expr_base.networkConstruction.TOM.RData")
geneTree = flashClust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
		deepSplit = 2, pamRespectsDendro = FALSE,
		minClusterSize = 20);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
\# Plot the dendrogram and the module colors underneath
pdf("$outDir/$expr_base.DendroAndColors.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
		dendroLabels = FALSE, hang = 0.03,
		addGuide = TRUE, guideHang = 0.05,
		main = "Gene dendrogram and module colors")
dev.off()

if(identical("$merge","Y")){
	merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.2,verbos=3)
	pdf("$outDir/$expr_base.DendroAndColors.after_merge.pdf")
	plotDendroAndColors(geneTree, cbind(dynamicColors,merged\$colors),
			c("Dynamic Tree Cut","Merged Dynamic"),dendroLabels = FALSE, hang = 0.03,
			addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
	dev.off()
	dynamicColors <- merged\$colors
}

save(datExpr,dynamicColors,adjacency,TOM,geneTree,file="$outDir/$expr_base.for_net_trait.RData")

save(datExpr,expr,dynamicMods,dynamicColors,geneTree,
			file = "$outDir/$expr_base.networkConstruction.RData")
q('no')
#;
	open O,">","$outDir/$expr_base.networkConstruction.R";
	print O $rscpt;
	close O;

	open O,">","$outDir/step2_net.sh";
	print O "R CMD BATCH $outDir/$expr_base.networkConstruction.R $outDir/$expr_base.networkConstruction.Rout\n";
	close O;

	`R CMD BATCH $outDir/$expr_base.networkConstruction.R $outDir/$expr_base.networkConstruction.Rout`;
}


