#!/usr/bin/perl -w
use strict;
my $genome_seq = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_internal_seq";

my %seq_hash = %{get_seq($genome_seq)};
my $input = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";
my $input1 = "$input"."_auto_mule";
my $output2 = "$input"."_Glimmer";
open(OUTPUT2,">$output2");
open(INPUT,$input);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $gene_id = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";
open(OUTPUT1,">temp_seq_japo");
print OUTPUT1 ">$gene_id\n$seq_hash{$gene_id}\n";
`rm MULE_gene/data/japo_MULE_Glimmer`;
`/wsu/home/fe/fe49/fe4960/GlimmerHMM/sources/glimmerhmm temp_seq_japo -d GlimmerHMM/trained_dir/rice -o MULE_gene/data/japo_MULE_Glimmer -g`;
open(INPUT1,"MULE_gene/data/japo_MULE_Glimmer");
<INPUT1>;
<INPUT1>;
while(my $line1=<INPUT1>){
chomp $line1;
print OUTPUT2 "$line1\n";
}
 
}

sub get_hash{
my $genelist = shift(@_);
my %hash;
open(INPUTL,$genelist);
while(my $lineL = <INPUTL>){
chomp $lineL;
my @info = split(/\s+/,$lineL);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $gene_id = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";
$hash{$gene_id}=1;
}
return(\%hash);
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
