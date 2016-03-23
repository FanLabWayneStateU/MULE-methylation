#!/usr/bin/perl -w
my @spe = ("japo","indi","rufi","niva","glab","bart","glum","meri","punc","brac","perr");
#my @spe = ("perr");
for (my $i=0; $i<=10; $i++){
my $goal_spe = $spe[$i];
my %hash;

for(my $j=0; $j<=$i-1; $j++){
my $target_spe = $spe[$j];
my $file = "MULE_gene/data/pack_MULE_candidate_$goal_spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_flank_list_$target_spe"."_ortho";

open(INPUT,"$file");

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$hash{$target_spe}->{$info[0]}=$info[1];
}
}


for(my $j=$i+1; $j<=$#spe; $j++){
my $target_spe = $spe[$j];
my $file = "MULE_gene/data/pack_MULE_candidate_$goal_spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_flank_list_$target_spe"."_ortho";

open(INPUT,"$file");

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$hash{$target_spe}->{$info[0]}=$info[1];
}
}

my $file_list = "MULE_gene/data/pack_MULE_candidate_$goal_spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_flank_list";

my $output = "MULE_gene/data/pack_MULE_candidate_$goal_spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho";

open(OUTPUT,">$output");
open(INPUT1,"$file_list");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
print OUTPUT "$info1[0]\t";
for(my $j=0; $j<=$i-1; $j++){
print OUTPUT "$spe[$j]:$hash{$spe[$j]}->{$info1[0]}\t";
}
for(my $j=$i+1; $j<=$#spe; $j++){
print OUTPUT "$spe[$j]:$hash{$spe[$j]}->{$info1[0]}\t";
}
print OUTPUT "\n";
} 

}

