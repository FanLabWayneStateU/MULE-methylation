#!/usr/bin/perl -w
use strict;
my $random_coor1 = "MULE_gene/data/japo_internal_MULE_c_coor_gene_no_te_maker_new1_mid200_deep_random_gene_nonTE";
#my $random_coor1 = "MULE_gene/data/niva_internal_MULE_c_coor_gene_no_te_maker_new1_mid200_deep_random_gene_nonTE";


open(OUTPUT1,">$random_coor1");
#my $genome_seq  = "MULE_gene/data/O.nivara_v1.0.pseudo.fasta";

my $genome_seq  = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};
my $genelist = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_non_TE_deep_no_transposase_non_TE_internal";
#my $genelist = "MULE_gene/data/maker_gff1/oryza_nivara.maker.gff_non_TE_deep_no_transposase_non_TE_internal";


open(INPUT,"$genelist");
while(my $line = <INPUT>){

chomp $line;

my @info = split(/\s+/,$line);
my $chr = substr($info[0],8);
if($chr =~ /^0(\d+)/){
$chr = "Chr"."$1";
}else{
$chr = "Chr$chr";
}
#my $chr = substr($info[0],5);

#if(abs($info[4]-$info[3])<=199){
#if(abs($info[4]-$info[3])<=164){
if(abs($info[4]-$info[3])<=160){

next;
}
my $real_name = substr($info[8],3,13);
my $length=abs($info[4]-$info[3])+1;
line1: my $random_number = int(rand($length));
#if(($info[3]+$random_number+199)>=$info[4]){
#if(($info[3]+$random_number+164)>=$info[4]){
if(($info[3]+$random_number+160)>=$info[4]){

goto line1;
}
$info[3] = $info[3]+$random_number;
#$info[4] = $info[3]+199;
#$info[4] = $info[3]+164;
$info[4] = $info[3]+160;

my ($q_old_start,$q_old_end,$q_old_chr) = ($info[3],$info[4],$chr);
open(OUTPUT3,">scaffold_old_no_te1");
print OUTPUT3 ">$q_old_chr\n$seq_hash{$q_old_chr}\n";

my $q_old_str="$info[6]";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_old_no_te1 $q_old_start $q_old_end $q_old_str > q_old_seq0_no_te1`;
open(INPUT1,"q_old_seq0_no_te1");
my $line1 = <INPUT1>;
chomp $line1;
my @info1 = split(//,$line1);
my $c_up;
if($q_old_str eq "+"){
	$c_up = $q_old_start;

	for(my $i=0; $i<=$#info1; $i++){
		if(($info1[$i] eq "C") || ($info1[$i] eq "c")){
	  	my $letter=uc($info1[$i]);	
		print OUTPUT1 "$real_name\t$q_old_chr\t$c_up\t$q_old_str\ttest\t$letter\n";
		}
		$c_up++;
	}
}elsif($q_old_str eq "-"){
	$c_up = $q_old_end;
	for(my $i=0; $i<=$#info1; $i++){
	  if(($info1[$i] eq "C")||($info1[$i] eq "c")){
		my $letter = uc($info1[$i]);
	  	print OUTPUT1 "$real_name\t$q_old_chr\t$c_up\t$q_old_str\ttest\t$letter\n";
	  }
        $c_up--;
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
