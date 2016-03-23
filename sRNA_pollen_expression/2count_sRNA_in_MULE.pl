#!/usr/bin/perl -w
use strict;
my $tissue = "TCP";
my @id = ("uniq_bwa0","two_more_bwa0");

my @file = ("MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon_genic_mule_non_auto1");

for my $id(@id){
my $rna = "MULE_gene/data/sRNA_tophat_"."$tissue/$tissue"."_sort_sRNA_list_"."$id";

my $sort_rna = "MULE_gene/data/sRNA_tophat_"."$tissue/$tissue"."_sort_sRNA_list_"."$id"."_sort";

`sort -k 2,2 -k 3,3n -k 4,4n $rna > $sort_rna`;


my %rna;

for my $file (@file){
open(INPUT2,"$file");
my $outputnew = "$sort_rna"."_mule_clean";

open(OUTPUT1,">$outputnew");

while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @start = split(/:/,$info2[2]);
my @end = split(/:/,$info2[3]);
my ($chr,$start,$end) = ($info2[1],$start[1],$end[0]);
my $name = "$info2[0]:$info2[1]:$start:$end:$start[2]";


my $internal_srna = get_rna_count($name,$chr,$start,$end,$sort_rna);
print OUTPUT1 "$name\t$internal_srna\n";

}
}
}


sub get_rna_count{
my ($gene,$chr,$cds_start,$cds_end,$file) = (@_);
`awk '{if(\$2 == "$chr"){print ;}}' $file | sort -k 3,3n -k 4,4n > sort_coor_new10`; 
open(INPUT1,"sort_coor_new10");
my $exp_flag=0;
my $exp_str="";
my $mule_overlap=0;
my $exon_overlap=0;
my %rna_total;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_start = $info1[2];
my $rna_end = $info1[3];

my $rna_name = $info1[0];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
 
		$exp_str .= "$rna_name:";
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$exp_str .= "$rna_name:";

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		$exp_str .= "$rna_name:";

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			$exp_str .= "$rna_name:";

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}
	}

}
line1: return($exp_str);

}
