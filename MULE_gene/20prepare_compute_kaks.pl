#!/usr/bin/perl -w
use strict;
my @target_spe = ("bart","brac","glum","meri","rufi","perr","japo","niva","glab","punc","indi");
my $spe = $ARGV[0];
my $singleton = "MULE_gene/data/$spe"."_MULE_gene_all_gene_seq";
my %single_seq =%{get_seq($singleton)};
my $mule_seq = "MULE_gene/data/$spe"."_MULE_internal_seq";
my %mule_seq = %{get_seq($mule_seq)};
my $blat_best = "MULE_gene/data/$spe"."_MULE_all_MULE.psl_best_sort";
my $parent_seq = "MULE_gene/data/$spe"."_MULE_gene_all_para_seq.fa";
open(OUTPUT1,">$parent_seq");
my $mule_seq2 = "MULE_gene/data/$spe"."_MULE_gene_all_mule_cds.fa";
open(OUTPUT2,">$mule_seq2");
my $kaks_list = "MULE_gene/data/$spe"."_MULE_gene_all_cds_list";
open(OUTPUT,">$kaks_list");
open(INPUT,$blat_best);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $temp_seq = "MULE_gene/data/$spe"."_temp_single";
my $temp_seq1 =  "MULE_gene/data/$spe"."_temp_single1";
open(OUTPUT3,">$temp_seq");
open(OUTPUT4,">$temp_seq1");
my $single_seq1 = $single_seq{$info[0]};
my $mule_seq1 = $mule_seq{$info[4]};
my $gene_start = $info[2]+1;
my $gene_end = $info[3];
my $mule_start = $info[6]+1;
my $mule_end = $info[7];
my $mule_str = $info[11];
my @mule_id = split(/:/,$info[4]);
my $mule_start1 = $mule_id[2]+$mule_start -1;
my $mule_end1 = $mule_id[2]+$mule_end -1;
my $new_mule_id = "$mule_id[0]:$mule_id[1]:$mule_start1:$mule_end1:$mule_str";
print OUTPUT "$info[0]\t$new_mule_id\n";
print OUTPUT2 ">$info[0]\n$single_seq1\n";

print OUTPUT3 ">$new_mule_id\n$mule_seq1\n";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $temp_seq $mule_start $mule_end $mule_str > $temp_seq1`;
open(INPUTt1,"$temp_seq1");
my $linet1 = <INPUTt1>;
chomp $linet1;
print OUTPUT1 ">$new_mule_id\n$linet1\n";
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
