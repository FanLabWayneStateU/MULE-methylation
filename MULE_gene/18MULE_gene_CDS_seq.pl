#!/usr/bin/perl -w
use strict;
my @target_spe = ("bart","brac","glum","meri","rufi","perr","japo","niva","glab","punc","indi");
for (my $k=0; $k<=10; $k++){
my $spe = $target_spe[$k];
my $cds_seq = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq";
my %cds_seq =%{get_seq($cds_seq)};
my $singleton =  "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_mRNA";
my $output = "MULE_gene/data/$spe"."_MULE_gene_all"."_gene_seq";

open(OUTPUT,">$output");
open(INPUT,"$singleton");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @info1 = split(/;/,$info[8]);
my $gene_name = substr($info1[0],3);
print OUTPUT ">$gene_name\n$cds_seq{$gene_name}\n";
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
