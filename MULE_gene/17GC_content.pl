#!/usr/bin/perl -w
use strict;
my @spe = ("japo","glum","indi","niva","perr","punc","rufi","bart","glab","meri","brac");

for my $spe (@spe){

#my $seq_file = "MULE_gene/data/maker_gff1/$spe".".maker.gff_non_TE_deep_no_transposase_non_TE_internal_gene_seq";


#my $seq_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_gene_seq";
my $seq_file = "MULE_gene/data/$spe"."_MULE_parental_seq_genome_blastn_gene_new1_exon30_parental_parental_seq";


#my $seq_file  = "MULE_gene/data/pack_MULE_$spe"."_gene_mule_part_exon30";


my $output = $seq_file."_gc";
open(OUTPUT,">$output");
my %hash_seq = %{get_seq($seq_file)};
for my $key (sort keys %hash_seq){
my $seq = $hash_seq{$key};
my $gc_content = get_gc($seq);
print OUTPUT "$key\t$gc_content\n";
}
}
sub get_gc{
my $seq = shift(@_);
my $new_seq = uc($seq);
my $gc_count=0;
my $non_N=0;
my @info =split(//,$new_seq);
for(my $i=0; $i<=$#info; $i++){
if(($info[$i] eq "C") || ($info[$i] eq "G")){
$gc_count++;
}
if(($info[$i] eq "A") || ($info[$i] eq "T") || ($info[$i] eq "C") || ($info[$i] eq "G")){
$non_N++;
}
}
my $gc_content=0;
if($non_N >0){
$gc_content = $gc_count/$non_N;

if(($non_N/($#info+1))<=0.1){
$gc_content=0;
}
if($non_N <=20){
$gc_content=0;
}
}
return($gc_content);
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
