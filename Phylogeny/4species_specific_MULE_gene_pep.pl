#!/usr/bin/perl -w
use strict;
my @spe = ("niva","bart","perr","glum","japo","indi","glab","rufi","brac","punc","meri"); 
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq_pep_screen";

my %pro = %{get_seq($file)};  
my $spe_list = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage";
my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
open(INPUT2,$coor);
my %coor;
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @name = split(/:/,$info2[0]);
my $id = "$name[0]:$name[1]:$name[2]:$name[3]";
$coor{$info2[1]}=$id;
}
my %spe_list;
open(INPUT1,"$spe_list");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
$spe_list{$info[0]}=1;
}
my $output1 = $spe_list."_pep_seq";
open(OUTPUT1,">$output1");
my $output = $spe_list."_gene_list";

open(OUTPUT,">$output");
for my $key (keys %pro){
my @info = split(/:/,$key);
my $id = "$info[0]:$info[1]:$info[2]:$info[3]:$info[4]";
if(exists $spe_list{$coor{$id}}){
print OUTPUT1 ">$key\n";
print OUTPUT1 "$pro{$key}\n";
print OUTPUT "$key\n";
}
}
}




sub get_seq{
my $fa_file = shift(@_);
my %hash;
open(INPUTfa,$fa_file);
my $linefa = <INPUTfa>;
chomp $linefa;
my $id='';
my $seq='';
if($linefa =~ /^>(\S+)/){
$id = "$1";
}
while(my $linefa = <INPUTfa>){
chomp $linefa;
if($linefa =~ /^>(\S+)/){
$hash{$id} = $seq;
$id = "$1";
$seq ="";
}else{
$seq .= "$linefa";
}
}
$hash{$id} = $seq;
print "done!\n";
return(\%hash);
}
