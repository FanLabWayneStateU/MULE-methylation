#!/usr/bin/perl -w
use strict;
my $spe_fa = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq = %{get_seq($spe_fa)};

#my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon_gff_genic";
my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic";

 
my $seq_file = $file."_gene_seq";
open(OUTPUT1,">$seq_file");
open(INPUT2,"$file");
while(my $line2 = <INPUT2>){
if($line2 =~ /mRNA/){
chomp $line2;
my @info3 = split(/\s+/,$line2);
my @info1 = split(/;/,$info3[8]);
my $id;
my @info2 = split(/:/,$info1[0]);
$id = substr("$info2[0]:$info2[1]:$info2[2]:$info2[3]:$info2[4]:$info2[5]",3);
print "$id\n";


my $chr = $info3[0];
my $start = $info3[3];
my $end = $info3[4];
my $str = $info3[6];
my $chr_seq = $seq{$chr};
open(OUTPUT, ">temp_chr1_new");
print OUTPUT ">$chr\n";
print OUTPUT "$chr_seq";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new $start $end $str > cds_seq1_new`;
open(INPUTt,"cds_seq1_new");
my $linet = <INPUTt>;
chomp $linet;

print OUTPUT1 ">$id\n";
print OUTPUT1 "$linet\n";
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
