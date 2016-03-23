#!/usr/bin/perl -w
use strict;
my @spe = ("bart","brac","glab","glum","perr","punc","rufi","indi","japo","niva","meri","long");
for my $spe (@spe){
my $file = "MULE_gene/data/auto_MULE_$spe"."_new_tir_internal_tblastn_new";
my $file2 = $file."_name";
`awk '{if(\$11 < 1e-9){split(\$2,a,":");if((a[4]-a[3]+1)>=3000){print \$2;}}}' $file | sort | uniq | cat > $file2`;
open(INPUT,$file2);
my %hash;
while(my $line = <INPUT>){
chomp $line;
$hash{$line}=1;

}
my $file4 = "MULE_gene/data/transposase_MULE_$spe"."_new_tir_internal_tblastn_new";
my $file5 = $file4."_auto_MULE";
`awk '{if((\$1 ~ /MuDR/)&&(\$11 < 1e-9)){split(\$2,a,":");if((a[4]-a[3]+1)>=3000){print \$2;}}}' $file4 | sort | uniq | cat > $file5`;
open(INPUT5,"$file5");
while(my $line5 = <INPUT5>){
chomp $line5;
$hash{$line5}=1;
}

my $file1 = "MULE_gene/data/maker_MULE_$spe"."_new_tir_internal_tblastn_only_new";

my $file3 = $file1."_name";

`awk '{if(\$11 < 1e-9){print \$2;}}' $file1 | sort | uniq | cat > $file3`;
my %hash1;
open(INPUT2,"$file3");
while(my $line = <INPUT2>){
chomp $line;
$hash1{$line}=1;
}
my $no_trans = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

open(INPUT1,"$no_trans");
my $packfile = $no_trans."_pack_mule_only";

open(OUTPUT1,">$packfile");
my $autofile = "$no_trans"."_auto_mule";
open(OUTPUT2,">$autofile");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $name = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";
my $name1 ="$info[0]:$info[1]:$start[1]:$end[0]:$start[2]";
 
if((!(exists $hash{$name}))&&(exists $hash1{$name})){

print OUTPUT1 "$line1\n";
}
if(exists $hash{$name}){
print OUTPUT2 "$line1\n";
}
}
} 

