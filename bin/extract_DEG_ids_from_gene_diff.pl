use strict;
use warnings;
use Getopt::Long;

my $usage="\n$0
--diff    \t<str>\tthe differential expression file in the format similar to the output of cufflinks
--outDir  \t<str>\tthe output directory
--padj    \t<str>\tthe fdr cut off to select differentially expressed genes (default 0.05)
--pv      \t<str>\tthe pvalue cut off to select differentially expressed genes. if this pra

This scripts is used to parse the gene_expr.diff file, output the following results:
1. differentially expressed gene ids in the format Sample1-VS-Sample2.up, here the up
   is denote expression higher in Sample2
\n";

my($fold,$diff,$outdir,$padj,$pv);

GetOptions(
		"diff:s"=>\$diff,
		"outDir:s"=>\$outdir,
		"padj:f"=>\$padj,
		"pv:f"=>\$pv
		);

die $usage if ! defined $diff;

$fold=2;
$padj ||= 0.05;

`mkdir -p $outdir` if ! -e "$outdir";

open I,$diff;
## GeneID  SampA   SampB   logFC   logCPM  PValue  AdjustP
my(@info,%res);
while(<I>){
 chomp;
 @info=split /\t/;
 next if /GeneID/;
 next if ($info[3] eq "NA" || $info[5] eq "NA" || $info[6] eq "NA");
 if(defined $pv){
 	$res{$info[1]."-VS-".$info[2]}{"up"}{$info[0]}++ if ($info[3] > log($fold)/log(2) && $info[5] < $pv);
 	$res{$info[1]."-VS-".$info[2]}{"down"}{$info[0]}++ if ($info[3] < -log($fold)/log(2) && $info[5] < $pv);
 }else{
	 if($info[3] > log($fold)/log(2) && $info[6] < $padj){
		 $res{$info[1]."-VS-".$info[2]}{"up"}{$info[0]}++;
	 }
	 if($info[3] < -log($fold)/log(2) && $info[6] < $padj){
		 $res{$info[1]."-VS-".$info[2]}{"down"}{$info[0]}++;
	 }
 }
}
close I;

open O2,">","$outdir/degNum";
print O2 "Up\tDown\n";

foreach my $samp (sort keys %res){
 open O,">","$outdir/$samp.up";
 print O join("\n",keys %{$res{$samp}{"up"}}),"\n";
 close O;
 open O,">","$outdir/$samp.down";
 print O join("\n",keys %{$res{$samp}{"down"}}),"\n";
 close O;
 print O2 "$samp\t";
 print O2 scalar keys %{$res{$samp}{"up"}},"\t";
 print O2 scalar keys %{$res{$samp}{"down"}},"\n";
 
}
close O2;

