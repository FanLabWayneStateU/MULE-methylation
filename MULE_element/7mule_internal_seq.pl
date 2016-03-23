#!/usr/bin/perl -w

use strict;
my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new";
#my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";


open(INPUT,"$file");
my $output = $file."_internal_seq";
#my $output = $file."_internal_seq_only";

open(OUTPUT,">$output");
my $genome_seq = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
open(OUTPUT1,">scaffold_new_new1");
print OUTPUT1 ">$info[1]\n";
print OUTPUT1 "$seq_hash{$info[1]}\n";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new_new1 $start[0] $end[1] $start[2] > middle_mid1`;

#`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new_new1 $start[1] $end[0] $start[2] > middle_mid1`;

open(INPUT1,"middle_mid1");
my $line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
print OUTPUT ">$info[0]:$info[1]:$start[0]:$end[1]:$start[2]\n";

print OUTPUT "$line1\n";
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
