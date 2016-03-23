#!/usr/bin/perl -w
use strict;
my $blast_out = "MULE_gene/data/distance/japo_prob_genome_single";
my %blast;
open(INPUT,$blast_out);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$blast{$info[0]} = substr($info[1],0,13);
}

my $prob_id = "MULE_gene/data/distance/japo_prob_id_list";

`grep ">" MULE_gene/data/distance/rice_prob_chr_all_new | awk '{split(\$1,a,">");split(\$2,b,"(");split(b[2],c,")");print a[2],c[1]}' | cat > $prob_id`;
my %prob_id;
open (INPUT1,"$prob_id");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
$prob_id{$info1[1]} = $info1[0];
}

for(my $i=1; $i<=12; $i++){
my $id=$i;
if($i<10){
$id = "0".$i;
}
my $file = "MULE_gene/data/distance/chr$id"."_prob_id";
my $output = $file."_gene";
open(OUTPUT,">$output");
open(INPUT3,"$file");
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
if($info3[0] =~ /^\d+/){

print OUTPUT "$info3[0]\t$info3[1]\t";
for(my $j=2;$j<=$#info3; $j++){
print OUTPUT "$info3[$j]:$blast{$prob_id{$info3[$j]}}\t";
}
print OUTPUT "\n";
}elsif($info3[0] =~ /^F/){
print OUTPUT "$info3[0]\t$info3[1]\t$info3[2]\t$info3[3]\t";
for(my $j=4;$j<=$#info3; $j++){
print OUTPUT "$info3[$j]:$blast{$prob_id{$info3[$j]}}\t";
}
print OUTPUT "\n";

}
}

}






