#! user/bin/perl -W
use strict;
my @spe = ("japo");
for my $spe (@spe){
my $auto_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_auto_mule";
my $screen_file = "MULE_gene/data/$spe"."_internal_MULE_c_coor_gene_new_$spe"."_methy_w_methy3_1_stat_genic1";

open(INPUT,$auto_mule);
my %auto;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $id = "$info[0]:$info[1]:$start[0]:$end[1]";
#my $id = "$info[0]:$info[1]:$start[1]:$end[0]";

$auto{$id}=1;
}
my $lineage = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage";
sub get_hash{
my $file = shift(@_);
my %hash;
open(INPUTf,"$file");
while(my $linef = <INPUTf>){
chomp $linef;
my @infof = split(/\s+/,$linef);
$hash{$infof[0]}=1;
}
return(\%hash);
}
my %lineage = %{get_hash($lineage)};
my $African = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_Asian_all";

#my $African = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_African";
my %African =  %{get_hash($African)};

my $AA = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_AA";
my %AA =  %{get_hash($AA)};
my $spe9 = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_BB_all";
my %spe9 =  %{get_hash($spe9)};


 

my $genic_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_mule";

my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
open(INPUT3,"$coor");
my %coor;
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
#$coor{substr($info3[0],0,length($info3[0])-2)}=substr($info3[1],0,length($info3[1])-2);
$coor{substr($info3[1],0,length($info3[1])-2)}=substr($info3[0],0,length($info3[0])-2);

}
my %genic_mule;
open(INPUT2,"$genic_mule");
while(my $line2 = <INPUT2>){
chomp $line2;
my @info = split(/:/,$line2);
my $id = "$info[0]:$info[1]:$info[2]:$info[3]";
#$genic_mule{$line2}=1;
#$genic_mule{$coor{$id}}=1;
$genic_mule{$id}=1;

}
my $output = $screen_file."_genic_non_auto_lineage_new";

open(OUTPUT,">$output");
my $output1 = $screen_file."_nongenic_non_auto_lineage_new";

open(OUTPUT1,">$output1");
my $output2 = $screen_file."_genic_non_auto_Asian_all_new";

open(OUTPUT2,">$output2");
my $output3 = $screen_file."_nongenic_non_auto_Asian_all_new";


open(OUTPUT3,">$output3");
my $output4 = $screen_file."_genic_non_auto_AA_new";

open(OUTPUT4,">$output4");
my $output5 = $screen_file."_nongenic_non_auto_AA_new";

open(OUTPUT5,">$output5");

my $output6 = $screen_file."_genic_non_auto_BB_all_new";

open(OUTPUT6,">$output6");
my $output7 = $screen_file."_nongenic_non_auto_BB_all_new";

open(OUTPUT7,">$output7");

open(INPUT1,$screen_file);
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
#my @start = split(/:/,$info1[2]);
#my @end = split(/:/,$info1[3]);
#my $id = "$info1[0]:$info1[1]:$start[0]:$end[1]:+";
#my $id = "$info1[0]:$info1[1]:$start[0]:$end[1]";
my $id = $info1[0];
#my $id1 = "$info1[0]:$info1[1]:$start[1]:$end[0]";
#my $id1 = $info1[0];
#my $id = $info1[0];
#my @id = split(/:/,$info1[1]);
my @id = split(/:/,$info1[0]);

my $id1 = "$id[0]:$id[1]:$id[2]:$id[3]";
if(exists $auto{$id1}){

#if(exists $auto{$id1}){
;
}else{
if($genic_mule{$id1}){
#if($genic_mule{$coor{$id1}}){
#if($genic_mule{$info1[0]}){
if(exists $lineage{$coor{$id1}}){
#if(exists $lineage{$id1}){

print OUTPUT "$line1\n";
}elsif(exists $African{$coor{$id1}}){
#}elsif(exists $African{$id1}){

print OUTPUT2 "$line1\n";
}elsif(exists $AA{$coor{$id1}}){
#}elsif(exists $AA{$id1}){

print OUTPUT4 "$line1\n";
}elsif(exists $spe9{$coor{$id1}}){
#}elsif(exists $spe9{$id1}){

print OUTPUT6 "$line1\n";
}
}else{
if(exists $lineage{$coor{$id1}}){
#if(exists $lineage{$id1}){

print OUTPUT1 "$line1\n";
}elsif(exists $African{$coor{$id1}}){
#}elsif(exists $African{$id1}){

print OUTPUT3 "$line1\n";
}elsif(exists $AA{$coor{$id1}}){
#}elsif(exists $AA{$id1}){

print OUTPUT5 "$line1\n";
}elsif(exists $spe9{$coor{$id1}}){
#}elsif(exists $spe9{$id1}){

print OUTPUT7 "$line1\n";
}
}
}
}

}
