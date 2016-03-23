#!/usr/bin/perl -w
use strict;
my $random_coor1 = "MULE_gene/data/japo_internal_MULE_c_coor_single_tir";
open(OUTPUT1,">$random_coor1");

my $genome_seq  = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};
#my $genelist = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final4";
my $genelist = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

open(INPUT,"$genelist");
while(my $line = <INPUT>){
chomp $line;

my @info = split(/\s+/,$line);
my $chr = $info[1];
my @name1 = split(/:/,$info[2]);
my @name2 = split(/:/,$info[3]);
my $real_name = "$info[0]:$info[1]:$name1[1]:$name2[0]:$name1[2]";
my ($q_old_start,$q_old_end,$q_old_chr) = ($name1[0],$name1[1],$chr);
open(OUTPUT3,">scaffold_old_tir");
print OUTPUT3 ">$q_old_chr\n$seq_hash{$q_old_chr}\n";

my $q_old_str="$name1[2]";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_old_tir $q_old_start $q_old_end $q_old_str > q_old_seq_tir0`;
open(INPUT1,"q_old_seq_tir0");
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


#########

my ($q_old_start1,$q_old_end1,$q_old_chr1) = ($name2[0],$name2[1],$chr);
open(OUTPUT3,">scaffold_old_tir1");
print OUTPUT3 ">$q_old_chr1\n$seq_hash{$q_old_chr1}\n";

my $q_old_str1="$name1[2]";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_old_tir1 $q_old_start1 $q_old_end1 $q_old_str > q_old_seq_tir1`;
open(INPUT1,"q_old_seq_tir1");
my $line1 = <INPUT1>;
chomp $line1;
my @info1 = split(//,$line1);

my $c_up;
if($q_old_str1 eq "+"){
	$c_up = $q_old_start1;

	for(my $i=0; $i<=$#info1; $i++){
		if(($info1[$i] eq "C") || ($info1[$i] eq "c")){
	  	my $letter=uc($info1[$i]);	
		print OUTPUT1 "$real_name\t$q_old_chr1\t$c_up\t$q_old_str\ttest\t$letter\n";
		}
		$c_up++;
	}
}elsif($q_old_str1 eq "-"){
	$c_up = $q_old_end1;
	for(my $i=0; $i<=$#info1; $i++){
	  if(($info1[$i] eq "C")||($info1[$i] eq "c")){
		my $letter = uc($info1[$i]);
	  	print OUTPUT1 "$real_name\t$q_old_chr1\t$c_up\t$q_old_str\ttest\t$letter\n";
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
