#!/usr/bin/perl -w
use strict;
my $random_coor1 = "MULE_gene/data/japo_internal_MULE_c_coor_gene_new";
#my $random_coor1 = "MULE_gene/data/japo_internal_MULE_c_coor_promoter_new";
#my $random_coor1 = "MULE_gene/data/niva_internal_MULE_c_coor_gene_new";

open(OUTPUT1,">$random_coor1");
#my $genome_seq  = "MULE_gene/data/O.nivara_v1.0.pseudo.fasta";

my $genome_seq  = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};
my $genelist = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_mRNA";
#y $genelist = "MULE_gene/data/pack_MULE_candidate_niva_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic_mRNA";

open(INPUT,"$genelist");
while(my $line = <INPUT>){

chomp $line;

my @info = split(/\s+/,$line);
if($info[2] eq "mRNA"){
my $chr = $info[0];
my @name = split(/;/,$info[8]);
#my $real_name = substr($info[8],3);
my $real_name =substr($name[0],3);
my ($q_old_start,$q_old_end,$q_old_chr) = ($info[3],$info[4],$chr);
#my $q_old_start = $info[3]-1000;
#my $q_old_end = $info[4]+1000;
#my $q_old_chr = $chr;
open(OUTPUT3,">scaffold_old_gene");
print OUTPUT3 ">$q_old_chr\n$seq_hash{$q_old_chr}\n";

my $q_old_str="$info[6]";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_old_gene $q_old_start $q_old_end $q_old_str > q_old_seq0_gene`;
open(INPUT1,"q_old_seq0_gene");
my $line1 = <INPUT1>;
chomp $line1;
my @info1 = split(//,$line1);
my $c_up;
if($q_old_str eq "+"){
	$c_up = $q_old_start;
#	for(my $i=0; $i<=499; $i++){
#	for(my $i=0; $i<=999; $i++){

	for(my $i=0; $i<=$#info1; $i++){
		if(($info1[$i] eq "C") || ($info1[$i] eq "c")){
	  	my $letter=uc($info1[$i]);	
		print OUTPUT1 "$real_name\t$q_old_chr\t$c_up\t$q_old_str\ttest\t$letter\n";
		}
		$c_up++;
	}
}elsif($q_old_str eq "-"){
	$c_up = $q_old_end;
#	for(my $i=0; $i<=499; $i++){
#	for(my $i=0; $i<=999; $i++){

	for(my $i=0; $i<=$#info1; $i++){
	  if(($info1[$i] eq "C")||($info1[$i] eq "c")){
		my $letter = uc($info1[$i]);
	  	print OUTPUT1 "$real_name\t$q_old_chr\t$c_up\t$q_old_str\ttest\t$letter\n";
	  }
        $c_up--;
        } 
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

sub get_coor{
my ($my_line,$flag)=(@_);
open(INPUT3,"$my_line");
my $line = <INPUT3>;
chomp $line;
my @info = split(/\s+/,$line);
my $chr = $info[0];
my $str = $info[6];
my $start = $info[3];
my $end = $info[4];

return($start,$end,$str,$chr);
}
