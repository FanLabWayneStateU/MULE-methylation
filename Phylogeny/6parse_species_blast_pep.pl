#!/usr/bin/perl -w
use strict;
my @spe = ("glab","bart","indi","japo","rufi","niva","glum","meri","punc","brac","perr");
for (my $i=0; $i<=$#spe; $i++){
my $focus_spe =$spe[$i];
my $genome_spe;
my %ortho;
my $output = "MULE_gene/data/pack_MULE_species_specific_"."$focus_spe"."_genome_pep_ortho2";
open(OUTPUT,">$output");
my $gene_list = "MULE_gene/data/pack_MULE_candidate_$focus_spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage";
my $pep_len = "$gene_list"."_pep_seq";
my %list = %{get_list($pep_len)};

  for (my $j=0; $j<=$i-1; $j++){
  $genome_spe = $spe[$j];
  my $file = "MULE_gene/data/pack_MULE_species_specific_"."$focus_spe"."_$genome_spe"."_genome_pep2";
  my %hash = %{get_pep($file,$pep_len)};
	  for my $key (sort keys %hash){
	  $ortho{$key}->{$genome_spe} = $hash{$key};
   	  }
  }
   
  for (my $j=$i+1; $j<=10; $j++){
  $genome_spe = $spe[$j];
  my $file = "MULE_gene/data/pack_MULE_species_specific_"."$focus_spe"."_$genome_spe"."_genome_pep2";
  my %hash = %{get_pep($file,$pep_len)};
	  for my $key (sort keys %hash){
	  $ortho{$key}->{$genome_spe} = $hash{$key};
   	  }

  }
  for my $gene ( sort keys %list){

	print OUTPUT "$gene\t";
	  for (my $j=0; $j<=$i-1; $j++){
	  $genome_spe = $spe[$j];
	      if(exists $ortho{$gene}->{$genome_spe}){
	      print OUTPUT "$genome_spe:$ortho{$gene}->{$genome_spe}\t";
	      }else{
	      print OUTPUT "$genome_spe:0\t";
              }
          }

	  for (my $j=$i+1; $j<=10; $j++){
	  $genome_spe = $spe[$j];
	      if(exists $ortho{$gene}->{$genome_spe}){
	      print OUTPUT "$genome_spe:$ortho{$gene}->{$genome_spe}\t";
	      }else{
	      print OUTPUT "$genome_spe:0\t";
              }
          }
	      print OUTPUT "\n";
   }
}

sub get_list {
my $file = shift(@_);
open(INPUT,$file);
my %hash;
while(my $line = <INPUT>){
chomp $line;
if($line =~ /^>(\S+)/){
my $info = "$1";
#my @info = split(/\s+/,$line);
$hash{$info}=1;
#$hash{$info[0]}=1;
}
}
return(\%hash);
}

sub get_pep{
my ($file,$pep) = (@_);
my %pep_len = %{get_len($pep)};
open(INPUT1,"$file");
my %already;
my %hash;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $iden = $info1[2];
my $coverage = $info1[3]/$pep_len{$info1[0]};
if(!(exists $already{$info1[0]})){
#if(($iden >=70) &&($coverage >= 0.7)){
#if(($iden >=50) &&($coverage >= 0.7)){
if($info1[10]< 1e-10){
$hash{$info1[0]}=$info1[1];
$already{$info1[0]}=1;
}
}
}
return(\%hash);
}

sub get_len{
my $fa_file = shift(@_);
my %hash;
open(INPUTfa,$fa_file);
my $linefa = <INPUTfa>;
chomp $linefa;
my $id='';
my $seq='';
if($linefa =~ /^>(\S+)/){
$id = "$1";
}
while(my $linefa = <INPUTfa>){
chomp $linefa;
if($linefa =~ /^>(\S+)/){
$hash{$id} = length($seq);
$id = "$1";
$seq ="";
}else{
$seq .= "$linefa";
}
}
$hash{$id} = length($seq);
print "done!\n";
return(\%hash);
}

