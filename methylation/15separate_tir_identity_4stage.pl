#!/usr/bin/perl -w
use strict;
my @spe = ("rufi","niva","indi","glab","bart","japo");

for my $spe (@spe){
my $tir_iden = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_tir_mafft_result";
my $output1 = $tir_iden."_lineage_genic_new";

#my $output1 = $tir_iden."_lineage_genic_jiang1";

open(OUTPUT1,">$output1");
my $output2 = $tir_iden."_Asian_all_genic_new";

#my $output2 = $tir_iden."_Asian_genic_jiang1";


open(OUTPUT2,">$output2");
my $output3 = $tir_iden."_AA_genic_new";

#my $output3 = $tir_iden."_AA_genic_jiang1";

open(OUTPUT3,">$output3");
my $output4 = $tir_iden."_BB_all_genic_new";

#my $output4 = $tir_iden."_BB_all_genic_jiang1";

open(OUTPUT4,">$output4");
my $output5 = $tir_iden."_Asian_all_nongenic_new";
#my $output5 = $tir_iden."_Asian_nongenic_jiang1";

open(OUTPUT5,">$output5");
my $output6 = $tir_iden."_AA_nongenic_new";
#my $output6 = $tir_iden."_AA_nongenic_jiang1";

open(OUTPUT6,">$output6");
my $output7 = $tir_iden."_BB_all_nongenic_new";
#my $output7 = $tir_iden."_BB_all_nongenic_jiang1";

open(OUTPUT7,">$output7");
my $output8 = $tir_iden."_lineage_nongenic_jiang1";

#my $output8 = $tir_iden."_lineage_nongenic_new";
open(OUTPUT8,">$output8");
my $auto_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_auto_mule";

open(INPUT,$auto_mule);
my %auto;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
#my $id = "$info[0]:$info[1]:$start[0]:$end[1]";
my $id = "$info[0]:$info[1]:$start[1]:$end[0]";

$auto{$id}=1;
}
my $genic_mule = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_mule";

my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
#my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_jiang";

open(INPUT3,"$coor");
my %coor;
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
$coor{substr($info3[0],0,length($info3[0])-2)}=substr($info3[1],0,length($info3[1])-2);
#$coor{substr($info3[1],0,length($info3[1])-2)}=substr($info3[0],0,length($info3[0])-2);

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


my %lineage = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage")};
my %African = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_Asian_all")};
#my %African = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_Asian")};
#my %African = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_African")};

#my %African = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_African_all")};

my %AA = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_AA")};
#my %spe9 = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_9spe")};
my %spe9 = %{get_hash("MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_BB_all")};




open(INPUT1,"$tir_iden");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my @tir1 = split(/:/,$info1[2]);
my @tir2 = split(/:/,$info1[4]);
#if($tir1[0]>$tir2[0]){
if($tir1[0]<$tir2[0]){

$info1[5] = $tir1[0];
}else{
$info1[5] = $tir2[0];
}
if(exists $coor{$info1[0]}){
if(!(exists $auto{$info1[0]})){
#if(exists $pack_list{$info1[0]}){
#if(exists $lineage{$info1[0]}){

#print OUTPUT1 "$info1[0]\t$info1[5]\n";
#}els
if(exists $genic_mule{$coor{$info1[0]}}){
if(exists $lineage{$info1[0]}){
print OUTPUT1 "$info1[0]\t$info1[5]\n";

}elsif(exists $African{$info1[0]}){
print OUTPUT2 "$info1[0]\t$info1[5]\n";
}elsif(exists $AA{$info1[0]}){
print OUTPUT3 "$info1[0]\t$info1[5]\n";
}elsif(exists $spe9{$info1[0]}){
print OUTPUT4 "$info1[0]\t$info1[5]\n";
}
}else{
if(exists $lineage{$info1[0]}){
print OUTPUT8 "$info1[0]\t$info1[5]\n";

}elsif(exists $African{$info1[0]}){
print OUTPUT5 "$info1[0]\t$info1[5]\n";
}elsif(exists $AA{$info1[0]}){
print OUTPUT6 "$info1[0]\t$info1[5]\n";
}elsif(exists $spe9{$info1[0]}){
print OUTPUT7 "$info1[0]\t$info1[5]\n";
}

}
}
} 
}
}
sub get_hash{
my $file = shift(@_);
my %hash;
open(INPUT,$file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @name = split(/:/,$info[0]);
my $name="$name[0]:$name[1]:$name[2]:$name[3]";
$hash{$name}=1;
}
return(\%hash);
}


