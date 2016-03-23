#!/usr/bin/perl -w
use strict;
my @tag = ("uniq_bwa0","two_more_bwa0");
for my $tag (@tag){
my @tissue = ("TCP","BCP","UNM","leaf","root","panicle","seedling","shoot","Leaf1","Root1","Callus1","dcl1","dcl3","rdr2","wt");
my @file = ("_Asian_all_new","_BB_all_new","_AA_new");

my $dir = "MULE_gene/data/sRNA";
for my $tmp (@file){
my $file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp";

my $new_file = "$file"."_sRNA_all_$tag"."_mutant_mule_mutant_new";


open(OUTPUT,">$new_file");
my %tissue;
for my $tissue(@tissue){
opendir(DIR,"$dir");
my $tissue_dir = $tissue."_processed";
my @sRNA_file = grep {-f && ($_ =~ /$tissue_dir\.txt_list$/)} map {"$dir/$_"} readdir(DIR);

closedir(DIR);
my %tissue_list; 
open(INPUT0,$sRNA_file[0]);
while(my $line0 = <INPUT0>){
chomp $line0;
my @info0 = split(/\s+/,$line0);
$tissue_list{$info0[0]}->{len} = $info0[2];
$tissue_list{$info0[0]}->{copy}= $info0[1];
}
my %hash;
my $srna_file1 = "MULE_gene/data/sRNA_tophat_$tissue/$tissue"."_sort_sRNA_list_$tag"."_sort_mule_clean";
open(INPUTs1,"$srna_file1");
my %gene_info;
while(my $lines1 = <INPUTs1>){
chomp $lines1;
my @infos1 = split(/\s+/,$lines1);
my $gene_name = $infos1[0];
$gene_info{$gene_name} = $infos1[1];
}

open(INPUT,$file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @name = split(/:/,$info[0]);
my @srna_string1 = split(/:/,$gene_info{$info[0]});
my ($new_line,%count);
for(my $i=0; $i<=$#srna_string1; $i++){
$new_line .= "$srna_string1[$i]:";
$count{$tissue_list{$srna_string1[$i]}->{len}} += $tissue_list{$srna_string1[$i]}->{copy};
}

my $line_count = "";
for my $key (sort keys %count){
$line_count .= ":$key:$count{$key}";
}

if($new_line){
$tissue{$info[0]}->{$tissue}->{line}="$tissue:$new_line";
$tissue{$info[0]}->{$tissue}->{count} = "$tissue"."$line_count";

}else{
$tissue{$info[0]}->{$tissue}->{line}="$tissue:0";
$tissue{$info[0]}->{$tissue}->{count} = "$tissue:0";

}
}
}
open(INPUT1,$file);
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
print  OUTPUT "$line1\t";
for my $tissue1(@tissue){
print OUTPUT "$tissue{$info1[0]}->{$tissue1}->{line}\t$tissue{$info1[0]}->{$tissue1}->{count}\t";
}
print OUTPUT "\n";
}
}
}
sub get_srna{
my ($name,$rna,$rna1,$tissue_list) = (@_);
my %tissue_list = %{$tissue_list};
my %count;
my $string="";
`grep "$name" $rna > rna1_two`;

open(INPUT2,"rna1_two");
my $line2 = <INPUT2>;
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @sRNA = split(/:/,$info2[9]);
for(my $i=0; $i<=$#sRNA; $i++){
$count{$tissue_list{$sRNA[$i]}->{len}} += $tissue_list{$sRNA[$i]}->{copy};
$string .= "$sRNA[$i]:";
}


`grep "$name" $rna1 > rna2_two`;

open(INPUT3,"rna2_two");
my $line3 = <INPUT3>;
chomp $line3;
my @info3 = split(/\s+/,$line3);
my @sRNA1 = split(/:/,$info3[9]);
for(my $i=0; $i<=$#sRNA1; $i++){
$count{$tissue_list{$sRNA1[$i]}->{len}} += $tissue_list{$sRNA1[$i]}->{copy};
$string .= "$sRNA1[$i]:";
}




return($string,\%count);
}
