#!/usr/bin/perl -w
use strict;
my @spe = ("japo","bart","glum","indi","niva","perr","punc","rufi","glab","meri","brac");
for my $spe (@spe){
#my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon_gff_exon_genic_gene_seq_list";
my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_gene_seq_list";

my $mule_coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
my %mule_coor;
open(INPUT1,"$mule_coor");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
$mule_coor{$info1[1]}=$info1[0];
}
my $expression_out = "$gff_file"."_mule_par";
open(OUTPUT,">$expression_out");
open(INPUT,$gff_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
#my $mule_name = substr($info[0],0,(length($info[0])-2));
my @mule_name = split(/:/,$info[0]);
my $mule_name = "$mule_name[0]:$mule_name[1]:$mule_name[2]:$mule_name[3]:$mule_name[4]";
my @coor1 = split(/:/,$mule_coor{$mule_name});
my $chr = $coor1[1];
my $cds_start = $coor1[2];
my $cds_end = $coor1[3];
my $cds_str = $coor1[4];
my $genename=$info[1];

my ($exp_string,$exp_flag) = get_exp($chr,$cds_start,$cds_end,$genename,$spe);
print OUTPUT "$info[0]\t$info[1]\t$exp_string\n";

}

}


sub get_exp{

my ($chr,$cds_start,$cds_end,$genename,$spe)=(@_);
#my $expression_file = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";
my $expression_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon_gff_exon_genic";
#`awk '{if((\$1=="$chr_name")&&(\$3=="exon")){print }}' $expression_file | sort -k 4,4n -k 5,5n > gene_expression_glab4`;
`grep "$genename" $expression_file | awk '{if(\$3=="mRNA"){print }}' | sort -k 4,4n -k 5,5n > gene_expression_glab4`;
#exit;
open(INPUT1,"gene_expression_glab4");
my $exp_flag=0;
my $exp_str="";
my $mule_overlap=0;
my $exon_overlap=0;
my $rna_info;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_str = $info1[6];
my $rna_start = $info1[3];
my $rna_end = $info1[4];
my @rna_name=split(/;/,$info1[8]);
my $rna_name = $rna_name[0];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		$exon_overlap = ($rna_end - $cds_start+1)/($rna_end - $rna_start+1);
		$mule_overlap = ($rna_end - $cds_start+1)/($cds_end - $cds_start+1);
 
		$rna_info = "$chr:$cds_start:$rna_end:$rna_str";
		$exp_str .="$rna_info\t";
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$exon_overlap =($cds_end - $cds_start + 1)/($rna_end - $rna_start + 1);
		$mule_overlap =1;
		$rna_info = "$chr:$cds_start:$cds_end:$rna_str";
	
		$exp_str .="$rna_info\t";

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		$mule_overlap = ($rna_end - $rna_start+1)/($cds_end - $cds_start +1);
		$exon_overlap = 1;
		$rna_info = "$chr:$rna_start:$rna_end:$rna_str";
		$exp_str .="$rna_info\t";

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			$mule_overlap = ($cds_end - $rna_start+1)/($cds_end - $cds_start + 1);
			$exon_overlap = ($cds_end - $rna_start+1)/($rna_end - $rna_start + 1);
			$rna_info = "$chr:$rna_start:$cds_end:$rna_str";
			$exp_str .="$rna_info\t";

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}
	}

}
line1: return($exp_str,$exp_flag);

}

sub overlap{

my ($cds_start,$cds_end,$rna_start,$rna_end)=(@_);
my $exp_flag=0;
my $exp_str="";
my $overlap;
my $rna_info;
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line2;

		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$overlap =1;
		goto line2;
		}elsif($rna_end <= $cds_start){
		$overlap=0;
		$exp_flag=0;
		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line2;


		}elsif(($rna_end >= $cds_end) &&($rna_start <= $cds_end)){
			$overlap = ($cds_end - $rna_start)/($cds_end - $cds_start + 1);
			$exp_flag=1;
			goto line2;

		}elsif(($rna_end > $cds_end) && ($rna_start > $cds_end)){
			$exp_flag=0;
			
		}
	}

line2: return($exp_flag);
}


