#!/usr/bin/perl -w
use strict;
my $spe_fa = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq = %{get_seq($spe_fa)};
my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list";
my $output = $file."_flank_seq";
open(OUTPUT,">$output");
my $output2 = $file."_flank_list";
open(OUTPUT1, ">$output2");
my $agi_annotation = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";
my $agi_annotation_mRNA = $agi_annotation."_gene";
my %gene;
my $Glimmer = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS";

my $Glimmer_mRNA = $Glimmer."_mRNA";
`awk '{if(\$3 == "mRNA"){print}}' $Glimmer > $Glimmer_mRNA`;
open(INPUT1,"$Glimmer_mRNA");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 =split(/\s+/,$line1);
my @name1 = split(/;/,$info1[8]);
my $gene_name1 = substr($name1[0],3);
$gene{$gene_name1}->{str} = $info1[6];
$gene{$gene_name1}->{start} = $info1[3];
$gene{$gene_name1}->{end} = $info1[4];
}

open(INPUT2,"$file");
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @info3 = split(/:/,$info2[0]);
my $chr = $info3[1];
my ($str,$start,$end,$pro_start_flank_start1,$pro_start_flank_end1,$pro_start_flank_start2,$pro_start_flank_end2,$pro_end_flank_start1,$pro_end_flank_end1,$pro_end_flank_start2,$pro_end_flank_end2);
open(OUTPUTt,">temp_chr1_new1");
print OUTPUTt ">$chr\n";
print OUTPUTt "$seq{$chr}\n";

if($info2[0] =~ /gene/){
my $gene_name = $info2[0];
$str = $gene{$gene_name}->{str};
$start = $gene{$gene_name}->{start};
$end = $gene{$gene_name}->{end};

}elsif($info2[0] =~ /Osjap/){
my $gene_name_new = $info2[0];
$str = $gene{$gene_name_new}->{str};
$start = $gene{$gene_name_new}->{start};
$end = $gene{$gene_name_new}->{end};

}

$pro_start_flank_start1 = $start -1000;
$pro_start_flank_end1 = $start -1;
$pro_start_flank_start2 = $start - 2000;
$pro_start_flank_end2 = $start - 1001;

$pro_end_flank_start1 = $end +1;
$pro_end_flank_end1 = $end + 1000;
$pro_end_flank_start2 = $end + 1001;
$pro_end_flank_end2 = $end + 2000;
my $seq_string="";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new1 $start $end $str > cds_seq1_new1`;
open(INPUTt,"cds_seq1_new1");
my $linet = <INPUTt>;
chomp $linet;
$linet = uc($linet);
print OUTPUT ">$info2[0]\n";
print OUTPUT "$linet\n";
print OUTPUT1 "$info2[0]\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new1 $pro_start_flank_start1 $pro_start_flank_end1 $str > cds_seq1_new1`;
open(INPUTt,"cds_seq1_new1");
my $linet = <INPUTt>;
chomp $linet;
$linet = uc($linet);
my $gene_id = "$info3[0]:$info3[1]:$pro_start_flank_start1:$pro_start_flank_end1";
print OUTPUT ">$gene_id\n";
print OUTPUT "$linet\n";
print OUTPUT1 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new1 $pro_start_flank_start2 $pro_start_flank_end2 $str > cds_seq1_new1`;
open(INPUTt,"cds_seq1_new1");
my $linet = <INPUTt>;
chomp $linet;
$linet = uc($linet);
my $gene_id = "$info3[0]:$info3[1]:$pro_start_flank_start2:$pro_start_flank_end2";
print OUTPUT ">$gene_id\n";
print OUTPUT "$linet\n";
print OUTPUT1 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new1 $pro_end_flank_start1 $pro_end_flank_end1 $str > cds_seq1_new1`;
open(INPUTt,"cds_seq1_new1");
my $linet = <INPUTt>;
chomp $linet;
$linet = uc($linet);
my $gene_id = "$info3[0]:$info3[1]:$pro_end_flank_start1:$pro_end_flank_end1";
print OUTPUT ">$gene_id\n";
print OUTPUT "$linet\n";
print OUTPUT1 "$gene_id\t";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new1 $pro_end_flank_start2 $pro_end_flank_end2 $str > cds_seq1_new1`;
open(INPUTt,"cds_seq1_new1");
my $linet = <INPUTt>;
chomp $linet;
$linet = uc($linet);
my $gene_id = "$info3[0]:$info3[1]:$pro_end_flank_start2:$pro_end_flank_end2";
print OUTPUT ">$gene_id\n";
print OUTPUT "$linet\n";
print OUTPUT1 "$gene_id\t";
print OUTPUT1 "\n";


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
