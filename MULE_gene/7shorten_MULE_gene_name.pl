#!/usr/bin/perl -W
use strict;
my @spe = ("indi","niva","rufi","bart","glab","meri","glum","perr","punc","brac","japo");
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_gene_seq";
my $output1 = $file."_list";
open(OUTPUT1,">$output1");
my $output = $file."_short";
open(OUTPUT,">$output");
open(INPUT,$file);
my ($name,$pre_name,$k);
while(my $line = <INPUT>){
chomp $line;
if($line =~ /^>(\S+)/){
my $long_name = $1;
my @info = split(/:/,"$long_name");
$name ="$info[0]:$info[1]:$info[2]:$info[3]:$info[4]";

if($pre_name eq $name){
$k++;
print OUTPUT ">$name:$k\n";
print OUTPUT1 "$name:$k\t$long_name\n";
}else{
$k=0;
print OUTPUT ">$name:$k\n";
print OUTPUT1 "$name:$k\t$long_name\n";
}
$pre_name = $name;
}else{
print OUTPUT "$line\n";
}
}
}
