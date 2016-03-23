#!/usr/bin/perl -w
use strict;
my @spe = ("bart","brac","glab","glum","meri","punc","rufi","indi","niva","japo","perr");
for my $spe (@spe){
my $seq_file = "MULE_gene/data/potential_MULE_repeatscout_$spe"."_20_jiang";
my $file = "MULE_gene/data/pack_MULE_candidate_"."$spe"."_jiang_repeatscout_new_tir_final3_name_5";
my $output = $seq_file."_de_nov_candidate";
my $output1 = "miss_MULE";
open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
open(INPUTa,$file);
my %mule_agi;
while(my $linea = <INPUTa>){
chomp $linea;
if(($linea =~ /R\=/) && (!($linea =~ /X/)) && (!($linea =~ /\?/))){
my @info = split(/\:/,$linea);
$mule_agi{$info[1]}=1;
}
}
for my $key (sort keys %seq_hash){
print $key;
if((exists $mule_agi{$key} )){
print "yes\n";
print OUTPUT ">$key\n";
print OUTPUT "$seq_hash{$key}\n";
}
}
}

sub get_seq{
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
$hash{$id} = $seq;
$id = "$1";
$seq ="";
}else{
$seq .= "$linefa";
}
}
$hash{$id} = $seq;
print "done!\n";
return(\%hash);
}

