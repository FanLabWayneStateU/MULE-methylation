#!/usr/bin/perl -w
use strict;
my @tag = ("uniq_bwa0","two_more_bwa0");

for my $tag (@tag){
my @tissue = ("TCP","BCP","UNM","leaf","root","panicle","seedling","shoot","Leaf1","Root1","Callus1","dcl1","dcl3","rdr2","wt");
my @file = ("_Asian_all_new","_AA_new","_BB_all_new");

my $dir = "MULE_gene/data/sRNA";
for my $tmp (@file){
my $file = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_methy_stat_w_genic"."$tmp";

my $new_file = "$file"."_sRNA_all_$tag"."_mutant_tir_mutant_new";

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
my $srna_file1 = "MULE_gene/data/sRNA_tophat_$tissue/$tissue"."_sort_sRNA_list_$tag"."_sort_tir_clean";
open(INPUTs1,"$srna_file1");
my %gene_info;
while(my $lines1 = <INPUTs1>){
chomp $lines1;
my @infos1 = split(/\s+/,$lines1);
my $gene_name = "$infos1[0]";
$gene_info{$gene_name} = $infos1[1];

}

my %exists;
open(INPUT,$file);
print "$file\n";
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @name = split(/:/,$info[0]);
my $coor_file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon_genic_mule_non_auto1";

`grep "$name[0]" $coor_file | grep "$name[1]" | grep "$name[2]" | grep "$name[3]" > coor_line`;
open(INPUTc,"coor_line");
my $linec = <INPUTc>;
chomp $linec;
my @infoc = split(/\s+/,$linec);
my @start = split(/:/,$infoc[2]);
my @end = split(/:/,$infoc[3]);
my $name1 = "$infoc[0]:$infoc[1]:$start[0]:$start[1]:$start[2]";
my $name2 = "$infoc[0]:$infoc[1]:$end[0]:$end[1]:$end[2]";
my @mule_name = ("$name1","$name2");
for my $mule_name (@mule_name){
if(! (exists $exists{$mule_name})){
$exists{$mule_name}=1;
my @srna_string1 = split(/:/,$gene_info{$mule_name});
print "@srna_string1\n";
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
$tissue{$mule_name}->{$tissue}->{line}="$tissue:$new_line";
$tissue{$mule_name}->{$tissue}->{count} = "$tissue"."$line_count";

}else{
$tissue{$mule_name}->{$tissue}->{line}="$tissue:0";
$tissue{$mule_name}->{$tissue}->{count} = "$tissue:0";

}

}
}
}
}
open(INPUT1,$file);
my %exists1;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my @name_mule = split(/:/,$info1[0]);
my $coor_file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon_genic_mule_non_auto1";

`grep "$name_mule[0]" $coor_file | grep "$name_mule[1]" | grep "$name_mule[2]" | grep "$name_mule[3]" > coor_line`;
open(INPUTc,"coor_line");
my $linec = <INPUTc>;
chomp $linec;
my @infoc = split(/\s+/,$linec);
my @start = split(/:/,$infoc[2]);
my @end = split(/:/,$infoc[3]);
my $name1 = "$infoc[0]:$infoc[1]:$start[0]:$start[1]:$start[2]";
my $name2 = "$infoc[0]:$infoc[1]:$end[0]:$end[1]:$end[2]";
my @mule_name = ("$name1","$name2");
for my $name_mule(@mule_name){
if(!(exists $exists1{$name_mule})){
$exists1{$name_mule}=1;
print  OUTPUT "$line1\t";
for my $tissue1(@tissue){
print OUTPUT "$tissue{$name_mule}->{$tissue1}->{line}\t$tissue{$name_mule}->{$tissue1}->{count}\t";
}
print OUTPUT "\n";
}
}
}
}
}
sub get_srna{
my ($name,$rna,$rna1,$tissue_list) = (@_);
my %tissue_list = %{$tissue_list};
my %count;
my $string="";
`grep "$name" $rna > rna1_tir_two`;

open(INPUT2,"rna1_tir_two");
my $line2 = <INPUT2>;
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @sRNA = split(/:/,$info2[9]);
for(my $i=0; $i<=$#sRNA; $i++){
$count{$tissue_list{$sRNA[$i]}->{len}} += $tissue_list{$sRNA[$i]}->{copy};
$string .= "$sRNA[$i]:";
}


`grep "$name" $rna1 > rna2_tir_two`;

open(INPUT3,"rna2_tir_two");
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
