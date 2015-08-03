#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--go      \t<str>\tgene ontology information. No header information GenID<tb>GOid<tb>Evid
--genDes  \t<str>\tgene description. No header information GenID<tb>GenDescription
--genChr  \t<str>\tgene chromosome information. GenID<tb>Chromosome
--taxid   \t<str>\ttaxonomy ID from NCBI taxonomy database
--genus   \t<str>\tgenus
--species \t<str>\tspecies
--outDir  \t<str>\toutput directory

\n";

my($go,$outDir,$genDes,$genChr,$taxid,$genus,$species);

GetOptions(
		"go:s"=>\$go,
		"outDir:s"=>\$outDir,
		"genDes:s"=>\$genDes,
		"genChr:s"=>\$genChr,
		"taxid:s"=>\$taxid,
		"genus:s"=>\$genus,
		"species:s"=>\$species
);

die $usage if !defined $go;

`mkdir -p $outDir` if ! -e $outDir;

my $rscp=qq#
library(AnnotationForge)
options(stringsAsFactors=FALSE)
fChr <- read.table("$genChr",sep="\\t")[,1:2]
colnames(fChr) <- c("GID","CHROMOSOME")

fGO <- read.table("$go",sep="\\t",quote = "")
fGO <- data.frame(fGO[,1:2],rep("IEA",dim(fGO)[1]))
colnames(fGO) <- c("GID","GO","EVIDENCE")

fSym <- read.table("$genDes",sep="\\t",quote = "")
fSym <- data.frame(fSym[,1],fSym[,1],fSym[,2])
fSym <- fSym[which(!is.na(fSym[,3])),]
colnames(fSym) <- c("GID","SYMBOL","GENENAME")

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
		version="0.1",
		maintainer="Pengcheng Yang <yangpc\@mail.biols.ac.cn>",
		author="Pengcheng Yang <yangpc\@mail.biols.ac.cn>",
		outputDir = "$outDir",
		tax_id="$taxid",
		genus="$genus",
		species="$species",
		goTable="go")

packge.name <- paste("org.",substr("$genus",1,1),"$species",".eg.db",sep="")
pckg <- paste("$outDir/",package.name,sep="")
\#install the packages
install.packages(pckg,repos=NULL)

makeDBPackage("GO_DB",
		affy=FALSE,
		prefix="lmi_GO",
		fileName="$go",
		baseMapType="gb",
		outputDir = tmpout,
		version="1.0.0",
		manufacturer = "Affymetrix",
		chipName = "Human Cancer G110 Array",
		manufacturerUrl = "http://www.affymetrix.com")

hcg110_IDs = system.file("extdata",
		"hcg110_ID",
		package="AnnotationDbi")
tmpout = "./"
	makeDBPackage("HUMANCHIP_DB",
			affy=FALSE,
			prefix="hcg110",
			fileName=hcg110_IDs,
			baseMapType="gb",
			outputDir = tmpout,
			version="1.0.0",
			manufacturer = "Affymetrix",
			chipName = "Human Cancer G110 Array",
			manufacturerUrl = "http://www.affymetrix.com")

GO_DB = createSimpleBimap(
		"go",
		"GID",
		"GO",
		org.Lmigratoria.eg.db:::datacache,
		"GO_DB",
		"org.Lmigratoria.eg.db"
		)
#;

open O,">","$outDir/create_annotation_package_for_bioconductor.R";
print O $rscp;
close O;


