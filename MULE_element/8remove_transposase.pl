#!/usr/bin/perl -w
use strict;
my @spe = ("bart","brac","glab","glum","japo","meri","perr","punc","rufi","indi","niva","long");

for my $spe (@spe){
my $infile  ="MULE_gene/data/transposase_MULE_"."$spe"."_new_tir_internal_tblastn_new";

my $file = $infile."_e9";

`awk '{if((!(\$1 ~ /MuDR/))&&(\$11 < 1e-9)){print \$2;}}' $infile | sort | uniq | cat > $file`;

open(INPUT,"$file");
my %trans;
while(my $line = <INPUT>){
chomp $line;
$trans{$line}=1;
}
my $pro_file = "MULE_gene/data/pack_MULE_candidate_"."$spe"."_jiang_repeatscout_new_tir_final5_new";

my $output = $pro_file."_no_transposase_e9";

open(OUTPUT,">$output");
open(INPUT1,"$pro_file");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $name = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";
if(exists $trans{$name}){
;
}else{
print OUTPUT "$line1\n";
}
}

}
