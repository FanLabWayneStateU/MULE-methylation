#!/usr/bin/perl -w
use strict;
my @spe = ("bart","glab","glum","japo","meri","perr","punc","rufi","indi","niva","brac");

for my $spe (@spe){
my $ortho_seq = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq_pep_screen";
my $output = "$ortho_seq"."_exon_num";
open(OUTPUT,">$output");
my $gff = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_exon";
open(INPUT,"$ortho_seq");
my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
my ($exon_num);
while(my $line = <INPUT>){
chomp $line;
if($line =~ />(\S+)/){
$line = "$1";
my @info = split(/\s+/,$line);
my $mule = $info[0];
if($mule =~ /gene/){
`grep "$mule" $gff | grep "exon" | wc > mule_exon_num`;
open(INPUT3,"mule_exon_num");
my $line3 = <INPUT3>;
chomp $line3;
my @info3 = split(/\s+/,$line3);
$exon_num = $info3[1];
}else{
my @info5 = split(/:/,$mule);
my $gene_name = $info5[-1];
my $maker_gff;
my $maker_dir = "MULE_gene/data/maker_gff1/";
if($spe eq "bart"){
$maker_gff = $maker_dir."oryza_barthii.maker.gff";
}elsif($spe eq "glum"){
$maker_gff = $maker_dir."oryza_glumaepatula.maker.gff";
}elsif($spe eq "glab"){
$maker_gff = $maker_dir."oryza_glaberrima.maker.gff";
}elsif($spe eq "perr"){
$maker_gff = $maker_dir."leersia_perrieri.maker.gff";
}elsif($spe eq "meri"){
$maker_gff = $maker_dir."oryza_meridionalis.maker.gff";
}elsif($spe eq "rufi"){
$maker_gff = $maker_dir."oryza_rufipogon.maker.gff";
}elsif($spe eq "punc"){
$maker_gff = $maker_dir."oryza_punctata.maker.gff";
}elsif($spe eq "japo"){
$maker_gff = $maker_dir."oryza_sativa_japonica.maker.gff";
}elsif($spe eq "long"){
$maker_gff = $maker_dir."oryza_longistaminata.maker.gff";
}elsif($spe eq "indi"){
$maker_gff = $maker_dir."oryza_sativa_indica.maker.gff";
}elsif($spe eq "niva"){
$maker_gff = $maker_dir."oryza_nivara.maker.gff";
}elsif($spe eq "brac"){
$maker_gff = $maker_dir."oryza_brachyantha.maker.gff";
}

`grep "$gene_name" $maker_gff | grep "exon" | wc > gene_exon_num`;
open(INPUTg,"gene_exon_num");
my $lineg = <INPUTg>;
chomp $lineg;
my @info6 = split(/\s+/,$lineg);
$exon_num = $info6[1];
}

print OUTPUT "$line\t$exon_num\n";

}
}
}


