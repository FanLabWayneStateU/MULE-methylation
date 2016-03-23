#!/usr/bin/perl -w
use strict;
my $file_mRNA = "MULE_gene/data/japo_MULE_parental_seq_genome_blastn_gene_new1_exon30_parental";


my $file = $file_mRNA;
open(INPUT,$file_mRNA);
my $output_all = $file."_phy_dist_all_new";

open(OUTPUT_all,">$output_all");
print OUTPUT_all "set\tmap\tphys\n";
my $output_all_list = $file."_phy_dist_all_list_new";

open(OUTPUT_all_list,">$output_all_list");

my %already;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i = 1; $i<=$#info; $i++){

my @coor = split(/:/,$info[$i]);
my $gene_name="$info[0]:$info[$i]";

my $post = ($coor[1]+$coor[2])/2;
my $chr_name;
if($coor[0] =~ /Chr(\d+)/){
$chr_name = $1;
}
if($chr_name < 10){
$chr_name = "0$chr_name";
}

#if(!(exists $already{$gene_name})){
print OUTPUT_all "\"Ojap\"\t\"OsjapChr$chr_name\"\t$post\n";
print OUTPUT_all_list "\"Ojap\"\t\"OsjapChr$chr_name\"\t$post\t$gene_name\n";

$already{$gene_name}=1;
#}
}
}













