#!/usr/bin/perl -w
use strict;
my @spe = ("O.barthii_v1.3","O.sativa.japonica_IRGSP_MSU_v7","O.brachyantha.v1.4","O.glaberrima.v1","O.glumaepatula.v1.5","O.meridionalis.v1.3","O.punctata_v1.2","Oryza_rufipogon_W1943_pseudo_v1.1","O.sativa.indica_9311","Lperr_v1.4","O.nivara_v1.0.pseudo");

my @spe_flank = ("bart","japo","brac","glab","glum","meri","punc","rufi","indi","perr","niva");

for(my $j=5; $j<=5; $j++){

my $spe=$spe_flank[$j];
 my $flank_seq = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list_flank_list"; 

my $query_genome = "MULE_gene/data/$spe[$j]".".fasta";
for(my $i=0; $i<=$j-1; $i++){
my $blat_genome = "MULE_gene/data/$spe[$i]".".fasta";
my $spe_name;
if($blat_genome =~ /bart/){
$spe_name = "bart";
}elsif($blat_genome =~ /japo/){
$spe_name = "japo";
}elsif($blat_genome =~ /brac/){
$spe_name = "brac";
}elsif($blat_genome =~ /glab/){
$spe_name = "glab";
}elsif($blat_genome =~ /glum/){
$spe_name = "glum";
}elsif($blat_genome =~ /long/){
$spe_name = "long";
}elsif($blat_genome =~ /meri/){
$spe_name = "meri";
}elsif($blat_genome =~ /punc/){
$spe_name = "punc";
}elsif($blat_genome =~ /rufi/){
$spe_name = "rufi";
}elsif($blat_genome =~ /indi/){
$spe_name = "indi";
}elsif($blat_genome =~ /perr/){
$spe_name = "perr";
}elsif($blat_genome =~ /niva/){
$spe_name = "niva";
}
my $blat_out = "MULE_gene/data/$spe"."_species_specific_blat_$spe_name"."_genome2.psl";
`perl 7MULE_gene_ortholog.pl $spe $spe_flank[$i] $query_genome $blat_genome $blat_out $flank_seq`;
}

for(my $i=$j+1; $i <=$#spe; $i++){
my $blat_genome = "MULE_gene/data/$spe[$i]".".fasta";
my $spe_name;
if($blat_genome =~ /bart/){
$spe_name = "bart";
}elsif($blat_genome =~ /japo/){
$spe_name = "japo";
}elsif($blat_genome =~ /brac/){
$spe_name = "brac";
}elsif($blat_genome =~ /glab/){
$spe_name = "glab";
}elsif($blat_genome =~ /glum/){
$spe_name = "glum";
}elsif($blat_genome =~ /long/){
$spe_name = "long";
}elsif($blat_genome =~ /meri/){
$spe_name = "meri";
}elsif($blat_genome =~ /punc/){
$spe_name = "punc";
}elsif($blat_genome =~ /rufi/){
$spe_name = "rufi";
}elsif($blat_genome =~ /indi/){
$spe_name = "indi";
}elsif($blat_genome =~ /perr/){
$spe_name = "perr";
}elsif($blat_genome =~ /niva/){
$spe_name = "niva";
}
my $blat_out = "MULE_gene/data/$spe"."_species_specific_blat_$spe_name"."_genome2.psl";
`perl 7MULE_gene_ortholog.pl $spe $spe_flank[$i] $query_genome $blat_genome $blat_out $flank_seq`;

}
}
