#!/usr/bin/perl -w
use strict;
my @spe=("japo");
for my $spe(@spe){
my $genic_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon_genic_mule_non_auto1";
my $nongenic_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon_nongenic_mule_non_auto1";
 
open(INPUT,$genic_mule);
my %genic;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
#my $id = "$info[0]:$info[1]:$start[1]:$end[0]:$start[2]";
my $id = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";

$genic{$id}=1;
}

open(INPUT1,"$nongenic_mule");
my %nongenic;
while(my $line = <INPUT1>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
#my $id = "$info[0]:$info[1]:$start[1]:$end[0]:$start[2]";
my $id = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";

$nongenic{$id}=1;
}
my $mutation = "MULE_gene/data/$spe"."_internal_MULE_c_coor_gene_new_$spe"."_methy_w_methy3_1_stat";

my $output1 = $mutation."_genic1";
my $output2 = $mutation."_nongenic1";
open(OUTPUT,">$output1");
open(OUTPUT1,">$output2");
open(INPUT2,"$mutation");
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @id = split(/:/,$info2[0]);
my $id = "$id[0]:$id[1]:$id[2]:$id[3]:$id[4]";
if(exists $genic{$id}){

#if(exists $genic{$info2[0]}){
print OUTPUT "$line2\n";
}elsif(exists $nongenic{$id}){

#}elsif(exists $nongenic{$info2[0]}){
print OUTPUT1 "$line2\n";
}
}
}
