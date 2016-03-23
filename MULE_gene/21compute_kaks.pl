#!/usr/bin/perl -w
use strict;
#my @spe = ("japo","rufi","bart","glum","meri","punc","brac","perr","indi","niva","glab");
my @spe = ("japo");
for my $spe (@spe){
my $query_seq = "MULE_gene/data/$spe"."_MULE_gene_all_mule_cds.fa";
my $gff = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS";
my $hit_seq = "MULE_gene/data/$spe"."_MULE_gene_all_para_seq.fa";
my $prob_loc = "MULE_gene/data/$spe"."_MULE_all_kaks_problem";
my $detail = "MULE_gene/data/$spe"."_MULE_all_detail";

my $kaks = "MULE_gene/data/$spe"."_MULE_all_kaks";

my $rel = "MULE_gene/data/$spe"."_MULE_all_rel";

my $blat_out = "MULE_gene/data/$spe"."_MULE_all_blatout";

my $pair_list = "MULE_gene/data/$spe"."_MULE_gene_all_cds_list";

`perl 22modified_gKaKs.pl -output test_out1 -query_seq $query_seq -gff $gff -hit_seq $hit_seq -chrom all -blat bin/x86_64-redhat-linux-gnu/blat -blast2seq blast-2.2.20/bin/bl2seq -codeml paml4.7/bin/codeml -problem_loc $prob_loc -detail $detail -kaks_file $kaks -spe 5 -rel $rel -blat_out $blat_out -pair_list $pair_list`;
} 


