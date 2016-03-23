#!/usr/bin/perl -w
use strict;
my @spe = ("bart","japo","brac","glab","glum","meri","punc","rufi","indi","perr","niva");
my @spe_gtf = ("oryza_barthii","oryza_sativa_japonica","oryza_brachyantha","oryza_glaberrima","oryza_glumaepatula","oryza_meridionalis","oryza_punctata","oryza_rufipogon","oryza_sativa_indica","leersia_perrieri","oryza_nivara");
for(my $i=0; $i<=$#spe; $i++){
my $file = "MULE_gene/data/transposase_MULE_$spe[$i]"."_protein_new_blastp";
my $new_file = "MULE_gene/data/transposase_MULE_$spe[$i]"."_protein_new_blastp_1e-9";

`awk '{if(\$11 < 1e-9){print \$2}}' $file > $new_file `;

my $file1 = "MULE_gene/data/Mu_transposase_MULE_$spe[$i]"."_protein_new_blastp";
my $new_file1 = "MULE_gene/data/Mu_transposase_MULE_$spe[$i]"."_protein_new_blastp_1e-9";

`awk '{if(\$11 < 1e-9){print \$2}}' $file1 > $new_file1 `;


open(INPUT,$new_file);
my %hash;
while(my $line = <INPUT>){
chomp $line;
$line = substr($line,0,13);

$hash{$line}=1;
}

open(INPUT1,$new_file1);
while(my $line1 = <INPUT1>){
chomp $line1;
$line1 = substr($line1,0,13);

$hash{$line1}=1;
}
my $no_te = "MULE_gene/data/maker_gff1/$spe_gtf[$i]".".maker.gff_non_TE_deep";

my $output = $no_te."_no_transposase";
open(OUTPUT,">$output");
open(INPUT1,"$no_te");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $genename = substr($info1[8],3,13);
if(exists $hash{$genename}){
}else{
print OUTPUT "$line1\n";
}
}
}
