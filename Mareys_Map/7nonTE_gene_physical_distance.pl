#!/usr/bin/perl -w
use strict;
my $file = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_non_TE_deep_no_transposase_non_TE_internal";
my $gff_file = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_gene";

open(INPUT,$file);
my $output_all = $file."_phy_dist_all_new";


open(OUTPUT_all,">$output_all");
print OUTPUT_all "set\tmap\tphys\n";
my $output_all_list = $file."_phy_dist_all_list_new";

open(OUTPUT_all_list,">$output_all_list");

my %already;
while(my $line = <INPUT>){
chomp $line;
`grep "$line" $gff_file > gff_line`;
open(INPUT1,"gff_line");
$line = <INPUT1>;
my @info = split(/\s+/,$line);
my @coor = split(/[;|\.]/,$info[8]);
my $gene_name=substr($coor[0],3);

my $post = ($info[3]+$info[4])/2;
my $chr_id;
if($info[0] =~ /(Chr\d+)/){
$chr_id = $1;
}
my $chr_name;
$chr_name = $chr_id;



if(!(exists $already{$gene_name})){
print OUTPUT_all "\"Ojap\"\t\"Osjap$chr_name\"\t$post\n";
print OUTPUT_all_list "\"Ojap\"\t\"Osjap$chr_name\"\t$post\t$gene_name\n";

$already{$gene_name}=1;
}

}













