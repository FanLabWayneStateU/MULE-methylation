#!/usr/bin/perl -w
use strict;
my $genome_seq = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};
#my $parental_seq = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon_gff_exon_genic_gene_seq_list_mule_par";
my $parental_seq = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_gene_seq_list_mule_par";

open(INPUTp,"$parental_seq");
my %parent;
while(my $linep = <INPUTp>){
chomp $linep;
my @infop =split(/\s+/,$linep);
my $temp_name = $infop[0];
for(my $k=2; $k<=$#infop; $k++){

if(! (exists $parent{$temp_name}->{count})){
$parent{$temp_name}->{count}=0;
}
$parent{$temp_name}->{array}->[$parent{$temp_name}->{count}] = $infop[$k];
$parent{$temp_name}->{count}++;
}
}
my $gene_file = "MULE_gene/data/pack_MULE_japo_gene_mule_part_exon30";


my %par;
open(OUTPUT,">$gene_file");
for my $gene (sort keys %parent){
if(exists $parent{$gene}->{array}){
my @parent = @{$parent{$gene}->{array}};
my $whole_seq="";
my $genename = $gene;
for my $parent (@parent){
my @genename = split(/:/,$parent);
my $chr=$genename[0];
my $chr_seq = $seq_hash{$chr};
open(OUTPUT1,">genome_chr0");
print OUTPUT1 ">$chr\n";
print OUTPUT1 "$chr_seq\n";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl genome_chr0 $genename[1] $genename[2] $genename[3] > parental_seq0`;
open(INPUT1,"parental_seq0");
my $seq = <INPUT1>;
chomp $seq;
print OUTPUT ">$genename\n";
print OUTPUT "$seq\n";
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
