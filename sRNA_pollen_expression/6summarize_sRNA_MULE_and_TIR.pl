#!/usr/bin/perl -w
use strict;
my @tag = ("uniq_bwa0");
my $output = "MULE_gene/data/japo_4_pack_mule_11_tissue_sRNA_species_test_all_japo_methy_w_uniq";

open(OUTPUT,">$output");

my @file = ("_Asian_all_new","_AA_new","_BB_all_new");

my @nt = ("24","21");
for my $nt (@nt){
print OUTPUT "$nt\n";
print OUTPUT "\t\ttotal\tTCP\tBCP\tUNM\tleaf\troot\tpanicle\tseedling\tshoot\tLeaf1\tRoot1\tCallus1\tdcl1\tdcl3\trdr2\twt\n";

for my $tmp (@file){
my $file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp"."_sRNA_all_uniq_bwa0"."_mutant_mule_mutant_new";

my $par_num = get_par_num($file); 
print OUTPUT "$tmp\t$par_num\t";
for(my $i=8; $i<=36; $i=$i+2){
my %hash;
for my $tag (@tag){

my $file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp"."_sRNA_all_$tag"."_mutant_mule_mutant_new";

my $id = '$'."$i";
`awk '{if($id ~ /$nt/){print \$1}}' $file  > count1`;

open(INPUT1,"count1");
while(my $line1 = <INPUT1>){
chomp $line1;
$hash{$line1}=1;
}
}
my $count = keys %hash;
print OUTPUT "$count\t";
}
print OUTPUT "\n";

$file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp"."_CHH_sRNA_all_uniq_bwa0"."_mutant_tir_mutant_new";

my $par_num = get_par_num($file); 
print OUTPUT "$tmp"."_tir\t$par_num\t";
for(my $i=8; $i<=36; $i=$i+2){
my %hash;
my $tir_num=0;
for my $tag (@tag){
$file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp"."_CHH_sRNA_all_$tag"."_mutant_tir_mutant_new";

my $id = '$'."$i";
`awk '{if($id ~ /$nt/){print \$1}}' $file  | wc | cat > count1`;
open(INPUT1,"count1");
my $line1 = <INPUT1>;
my @info1 = split(/\s+/,$line1);
$tir_num += $info1[1];
}
print OUTPUT "$tir_num\t";
}
print OUTPUT "\n";



}
}


sub get_count{
my $file1 = shift(@_);
open(INPUT1,"$file1");
my $total_num=0;
if(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
$total_num =$info1[1];
}
return($total_num);
}

sub get_par_num{
my ($gene_list) = shift(@_);
open(INPUT2,"$gene_list");
my $count=0;
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
$count++;
}
return($count);
}


