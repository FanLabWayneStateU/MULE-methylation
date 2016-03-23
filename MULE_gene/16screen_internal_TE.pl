#!/usr/bin/perl -w
use strict;
my $gene_gff = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_non_TE_deep_no_transposase";

my $output = $gene_gff."_non_TE_internal";
open(OUTPUT,">$output");
open(INPUT,"$gene_gff");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $chr = substr($info[0],8);
if($chr =~ /^0/){
$chr = substr($chr,1);
}
my $start = $info[3];
my $end = $info[4];
$chr = "Chr$chr";
my $te = check_te($chr,$start,$end);
if($te ==0){
print OUTPUT "$line\n";
}
}

sub check_te{
my ($chr,$cds_start,$cds_end) = (@_);
my $te_gff = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";

`awk '{split(\$2,a,":");if(a[2]=="$chr"){print a[2],a[3],a[4],a[5]}}' $te_gff | sort -k 2,2n -k 3,3n > gene_te_japo`;
open(INPUT1,"gene_te_japo");
#exit;
my $exp_flag=0;
my $exp_str="";
my $mule_overlap=0;
my $exon_overlap=0;
my $rna_info;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_str = $info1[3];
my $rna_start = $info1[1]-1000;
my $rna_end = $info1[2]+1000;
#if($rna_str eq $cds_str){
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line1;
	#	$exon_overlap = ($rna_end - $cds_start+1)/($rna_end - $rna_start+1);
	#	$mule_overlap = ($rna_end - $cds_start+1)/($cds_end - $cds_start+1);
 
	#	$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$rna_end";
	#	$exp_str .="$rna_info\t";
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		goto line1;
	#	$exon_overlap =($cds_end - $cds_start + 1)/($rna_end - $rna_start + 1);
	#	$mule_overlap =1;
 	#	$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$rna_end";
	#	$exp_str .="$rna_info\t";

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line1;
	#	$mule_overlap = ($rna_end - $rna_start+1)/($cds_end - $cds_start +1);
	#	$exon_overlap = 1;
	#	$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$rna_end";
	#	$exp_str .="$rna_info\t";

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			goto line1;
	#		$mule_overlap = ($cds_end - $rna_start+1)/($cds_end - $cds_start + 1);
	#		$exon_overlap = ($cds_end - $rna_start+1)/($rna_end - $rna_start + 1);

	#		$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$rna_end";
	#		$exp_str .="$rna_info\t";

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}
	}
#}

}
line1: return($exp_flag);

}

