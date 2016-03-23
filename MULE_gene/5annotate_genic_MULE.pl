#!/usr/bin/perl -w
use strict;
my @spe = ("bart","brac","glab","glum","indi","japo","meri","niva","perr","punc","rufi");
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_"."$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq_pep_screen";

my %spe_list = %{get_list($file)};
sub get_list{
my $file = shift(@_);
open(INPUT,$file);
my %hash;
while(my $line = <INPUT>){
if($line =~ /^>(\S+)/){
my $id = "$1";
$hash{$id} = 1;
}
}
return(\%hash);

}
#my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS";
my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon";

my $output = $gff_file."_genic";
my $output1 = $gff_file."_genic_mule";


open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
open(INPUT1,$gff_file);
while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
my @info1 = split(/;/,$info[8]);
my $id;
my @info2;
if($info1[0] =~ /gene/){
@info2 = split(/:/,$info1[0]);
$id = substr("$info2[0]:$info2[1]:$info2[2]:$info2[3]:$info2[4]:$info2[5]",3); 
}else{
@info2 = split(/:/,$info1[0]);
$id = substr("$info2[0]:$info2[1]:$info2[2]:$info2[3]:$info2[4]:$info2[5]",3); 

}
if(exists $spe_list{$id}){
print OUTPUT "$line1\n";
print OUTPUT1 "$id\n";

}
}
}
