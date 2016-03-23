#!/usr/bin/perl -w
use strict;
my $gff_file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";

my $expression_out = "$gff_file"."_involved_exon30";

open(OUTPUT,">$expression_out");
open(INPUT,$gff_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @name = split(/:/,$info[1]);
my $cds_start = $name[2];
my $cds_end = $name[3];


my $cds_str = "+";
my $chr_num = substr($name[1],3);
my $chr_name;
if($chr_num < 10){
$chr_name = "OsjapChr0$chr_num";
}else{
$chr_name = "OsjapChr$chr_num";
}

my ($exp_string,$exp_flag) = get_exp($chr_name,$cds_start,$cds_end,$cds_str);
print OUTPUT "$info[0]\t$info[1]\t$exp_string\n";
}


sub get_exp{

my ($chr_name,$cds_start,$cds_end,$cds_str)=(@_);
my $expression_file = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";

`awk '{if((\$1=="$chr_name")&&(\$3=="exon")){print }}' $expression_file | sort -k 4,4n -k 5,5n > gene_expression_japo1`;

open(INPUT1,"gene_expression_japo1");
my $exp_flag=0;
my $exp_str="";
my $mule_overlap=0;
my $exon_overlap=0;
my $rna_info;
my %exon;
my %mule;
my %exon_str;
my %exon_len;
my %exon_strand;
my $mule_len = $cds_end - $cds_start + 1;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_str = $info1[6];
my $rna_start = $info1[3];
my $rna_end = $info1[4];
my @rna_name=split(/;/,$info1[8]);
my $rna_name = substr($rna_name[0],7,length($rna_name[0])-7);
$exon_len{$rna_name} += abs($info1[4] - $info1[3]) +1;
$exon_strand{$rna_name}="$info1[6]";

	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		$exon_overlap = ($rna_end - $cds_start+1)/($rna_end - $rna_start+1);
		$mule_overlap = ($rna_end - $cds_start+1)/($cds_end - $cds_start+1);
 		$exon{$rna_name} += $rna_end - $cds_start+1;
		$mule{$rna_name} += $rna_end - $cds_start+1;
		$exon_str{$rna_name} .= "$cds_start"."_"."$rna_end".":";
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$exon{$rna_name} += $cds_end - $cds_start + 1;
		$mule{$rna_name} += $cds_end - $cds_start + 1;
		$exon_str{$rna_name} .= "$cds_start"."_"."$cds_end".":";


		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		$exon{$rna_name} += $rna_end- $rna_start + 1;
		$mule{$rna_name} += $rna_end - $rna_start + 1;
		$exon_str{$rna_name} .= "$rna_start"."_"."$rna_end".":";


		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			$exon{$rna_name} += $cds_end - $rna_start + 1;
			$mule{$rna_name} += $cds_end - $rna_start + 1;
			$exon_str{$rna_name} .= "$rna_start"."_"."$cds_end".":";


		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
		;
		}
	}

}
line1: my $exon_overlap=0;
my $mule_overlap=0;

for my $key (sort keys %exon){
$exon_overlap = $exon{$key} /$exon_len{$key};
$mule_overlap = $mule{$key} / $mule_len;
#if($exon_overlap >= 0.1){
if($exon_overlap >= 0.3){

$rna_info = "$key:gene_"."$exon_overlap".":mule_"."$mule_overlap";
$exp_str .= "$rna_info:$exon_strand{$key}:$exon_str{$key}\t";
}
}
return($exp_str,$exp_flag);

}
