#!/usr/bin/perl -w
use strict;
my @spe = ("bart","brac","glab","glum","indi","japo","meri","niva","perr","punc","rufi");
for my $spe(@spe){
my $file = "MULE_gene/data/$spe"."_MULE_parental_seq_genome_blastn_gene_new1_exon30_parental_gene_structure";
my $file2 = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_gene_seq_list";
open(INPUT3,"$file2");
my %hash;
my %already;
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
if($info3[1] =~ /gene/){
$hash{$info3[0]}=$info3[1];
}else{
my @gene_name;

if($spe ne "long"){
@gene_name = split(/[:|\.]/,$info3[1]);
if(exists $already{$gene_name[-2]}){
;
}else{
$hash{$info3[0]} = $gene_name[-2];
$already{$gene_name[-2]}=1;
}

}elsif($spe eq "long"){
@gene_name = split(/[:|\-]/,$info3[1]);
$hash{$info3[0]} = $gene_name[-2];
if(exists $already{$gene_name[-2]}){
;
}else{
$hash{$info3[0]} = $gene_name[-2];
$already{$gene_name[-2]}=1;
}

}
}

my %spe;
open(INPUT,$file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i=2; $i<=$#info; $i++){
if(exists $hash{$info[0]}){
if(! (exists $spe{$hash{$info[0]}}->{count})){
$spe{$hash{$info[0]}}->{count}=0;
}
my @info2 = split(/[:|=]/,$info[$i]);
$spe{$hash{$info[0]}}->{gene}->{$info2[1]}=1;
$spe{$hash{$info[0]}}->{count}++;
}
}

}
#}


my $output = $file."_parental_gene_num";

open(OUTPUT,">$output");

	for my $key (sort keys %spe){
	my $num = keys (%{$spe{$key}->{gene}});
	print OUTPUT "$key\t$num\n";
	}
}
}









