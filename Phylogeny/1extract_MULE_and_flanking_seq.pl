#!/usr/bin/perl -w
use strict;
my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

my $output = $file."_flank_seq";
open(OUTPUT,">$output");
my $output2 = $file."_flank_list";
open(OUTPUT2,">$output2");
open(INPUT,$file);
my $genome_seq = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};

while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);

my $pro_start = $start[1];
my $pro_end = $end[0];
my $pro_start_flank_start1 = $pro_start-1000;
my $pro_start_flank_end1=$pro_start-1;
my $pro_start_flank_start2 = $pro_start-2000;
my $pro_start_flank_end2=$pro_start-1001;
my $pro_end_flank_start1 = $pro_end+1;
my $pro_end_flank_end1=$pro_end+1000;
my $pro_end_flank_start2 = $pro_end+1001;
my $pro_end_flank_end2=$pro_end+2000;


open(OUTPUT1,">scaffold_new11");
print OUTPUT1 ">$info[1]\n";
print OUTPUT1 "$seq_hash{$info[1]}\n";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new11 $pro_start $pro_end + > middle11`;
open(INPUT1,"middle11");
my $line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
my $gene_id = "$info[0]:$info[1]:$pro_start:$pro_end";
print OUTPUT ">$gene_id\n";
print OUTPUT "$line1\n";
print OUTPUT2 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new11 $pro_start_flank_start1 $pro_start_flank_end1 + > middle11`;
open(INPUT1,"middle11");
$line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
$gene_id="$info[0]:$info[1]:$pro_start_flank_start1:$pro_start_flank_end1";
print OUTPUT ">$gene_id\n";
print OUTPUT "$line1\n";
print OUTPUT2 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new11 $pro_start_flank_start2 $pro_start_flank_end2 + > middle11`;
open(INPUT1,"middle11");
$line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
$gene_id="$info[0]:$info[1]:$pro_start_flank_start2:$pro_start_flank_end2";
print OUTPUT ">$gene_id\n";
print OUTPUT "$line1\n";
print OUTPUT2 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new11 $pro_end_flank_start1 $pro_end_flank_end1 + > middle11`;
open(INPUT1,"middle11");
$line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
$gene_id="$info[0]:$info[1]:$pro_end_flank_start1:$pro_end_flank_end1";
print OUTPUT ">$gene_id\n";
print OUTPUT "$line1\n";
print OUTPUT2 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new11 $pro_end_flank_start2 $pro_end_flank_end2 + > middle11`;
open(INPUT1,"middle11");
$line1 = <INPUT1>;
chomp $line1;
$line1 = uc($line1);
$gene_id="$info[0]:$info[1]:$pro_end_flank_start2:$pro_end_flank_end2";
print OUTPUT ">$gene_id\n";
print OUTPUT "$line1\n";
print OUTPUT2 "$gene_id\n";


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
