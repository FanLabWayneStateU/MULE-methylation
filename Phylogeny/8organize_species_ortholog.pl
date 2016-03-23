#!/usr/bin/perl -w
use strict;
my @spe = ("bart","japo","indi","rufi","glab","niva","glum","meri","punc","brac","perr");
for(my $i=0; $i<=10; $i++){
my $spe1 = $spe[$i];
my %hash;

for (my $j=0; $j<=$i-1; $j++){
my $spe = $spe[$j];
my $file = "MULE_gene/data/pack_MULE_candidate_$spe1"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list_flank_list"."_$spe"."_ortho1";

open(INPUT,"$file");

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$hash{$spe}->{$info[0]}=$info[1];
}
}

for (my $j=$i+1; $j<=10; $j++){
my $spe = $spe[$j];

my $file = "MULE_gene/data/pack_MULE_candidate_$spe1"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list_flank_list"."_$spe"."_ortho1";


open(INPUT,"$file");

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$hash{$spe}->{$info[0]}=$info[1];
}
}


my $file_list = "MULE_gene/data/pack_MULE_candidate_$spe1"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list_flank_list";


my $output = "MULE_gene/data/pack_MULE_candidate_$spe1"."_jiang_repeatscout_final5_new_no_transposase_e9_11_spe_ortho_species";

open(OUTPUT,">$output");
open(INPUT1,"$file_list");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
print OUTPUT "$info1[0]\t";
for (my $j=0; $j <=$i-1; $j++){
my $spe1 = $spe[$j];
print OUTPUT "$spe1:$hash{$spe1}->{$info1[0]}\t";
}

for (my $j=$i+1; $j <=10; $j++){
my $spe1 = $spe[$j];
print OUTPUT "$spe1:$hash{$spe1}->{$info1[0]}\t";
}



print OUTPUT "\n";
} 

}
