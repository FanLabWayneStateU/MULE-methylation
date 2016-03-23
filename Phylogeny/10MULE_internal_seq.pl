#!/usr/bin/perl -w
use strict;
my @target_spe = ("bart","brac","glum","meri","rufi","perr","japo","niva","glab","punc","indi");
for (my $k=0; $k<=10; $k++){
my $spe = $target_spe[$k];
#my $cds_seq = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq";
#my %cds_seq =%{get_seq($cds_seq)};
#my $singleton = "MULE_gene/data/$spe"."_MULE_gene_singleton";
#my $output = $singleton."_gene_seq";
#open(OUTPUT,">$output");
#open(INPUT,"$singleton");
#while(my $line = <INPUT>){
#chomp $line;
#my @info = split(/\s+/,$line);
#print OUTPUT ">$info[0]\n$cds_seq{$info[0]}\n";
#}
my $mule_seq = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_internal_seq";
my %mule_seq = %{get_seq($mule_seq)};
my $all_gene_mRNA = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";

my $output1 = "MULE_gene/data/$spe"."_MULE_internal_seq";

open(OUTPUT1,">$output1");
open(INPUT1,$all_gene_mRNA);
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $id = $info1[1];
print OUTPUT1 ">$id\n$mule_seq{$id}\n";

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
