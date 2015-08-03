#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gdiff   \t<str>\tgene_diff file in the format: GeneID<tb>SampA<tb>SampB<tb>logFC<tb>logRPKM<tb>PValue<tb>FDR
--outDir  \t<str>\toutput directory

NOTE: the output file with the header information  GeneID<tb>Pvalue<tb>logFC<tb>SignedLogPv<tb>rank
SignedLogPv = -log(Pv)*logFC/abs(logFC)
the pvalue with 0 will be set to 1e-305 for logPv calculation
The rank is increased with the latter sample expressed higher

\n";

my($gdiff,$outDir);

GetOptions(
		"gdiff:s"=>\$gdiff,
		"outDir:s"=>\$outDir
);

die $usage if !defined $outDir;

open I,$gdiff;
my(@info,%gstat);
while(<I>){
	next if /^GeneID/;
	@info=split /\t/;
	my $pv=$info[5];$pv=1e-305 if $pv == 0;
	my($logpv);
	if($info[3]==0){
		$logpv=0;
	}else{
		$logpv = -log($pv) * $info[3]/abs($info[3]);
	}
	$logpv = 0 if $pv==1;
	$gstat{$info[1]}{$info[2]}{$info[0]}=$info[5]."\t".$info[3]."\t".$logpv;
}
close I;

`mkdir -p $outDir` if ! -e $outDir;
open O2,">","$outDir/gstat.list";
open O3,">","$outDir/gstat_rnk.list";
my($id1,$id2);
foreach $id1 (sort keys %gstat){
	foreach $id2 (sort keys %{$gstat{$id1}}){
		open O,">","$outDir/$id1"."_VS_$id2";
		map{print O "$_\t$gstat{$id1}{$id2}{$_}\n"}(sort keys %{$gstat{$id1}{$id2}});
		close O;
		print O2 $id1."_VS_$id2\t$outDir/$id1"."_VS_$id2\n";
		print O3 $id1."_VS_$id2\t$outDir/$id1"."_VS_$id2.rnk\n";
	}
}
close O2;


my $rscp=qq#
glist = read.table("$outDir/gstat.list",sep="\\t",stringsAsFactors=FALSE)
for(i in glist[,2]){
	a = read.table(i,stringsAsFactors=FALSE)
	a.res = data.frame(a,rank(a[,4],ties.method = "random"))
	i.rnk = paste(i,".rnk",sep="")
	write.table(a.res,file=i.rnk,sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	file.copy(from=i.rnk,to=i,overwrite = TRUE)
	file.remove(i.rnk)
	write.table(a.res[,c(1,5)],file=i.rnk,sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
#;

open O,">","$outDir/add_rank.R";
print O $rscp;
close O;
`R CMD BATCH $outDir/add_rank.R $outDir/add_rank.Rout`;

