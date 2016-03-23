#!/usr/bin/perl -w
use strict;
my $random_coor1 = "MULE_gene/data/japo_internal_MULE_c_coor_par_new_masked_gene_mule_parental_seq_new1";

open(OUTPUT1,">$random_coor1");

my $genome_seq  = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};


my $parlist = "MULE_gene/data/japo_MULE_parental_seq_genome_blastn_gene_new1_exon30_parental";

open(INPUT,"$parlist");
while(my $line = <INPUT>){

chomp $line;

my @info = split(/\s+/,$line);

for(my $i=1; $i<=$#info; $i++){
my @infos=split(/:/,$info[$i]);
my @chr = split(/:/,$info[$i]);
my $chr = $chr[0];
my $real_name = "$info[0]:$info[$i]";

my ($q_old_start,$q_old_end,$q_old_chr) = ($infos[1],$infos[2],$chr);
open(OUTPUT3,">scaffold_old_par0_1");
print OUTPUT3 ">$q_old_chr\n$seq_hash{$q_old_chr}\n";

my $q_old_str="$infos[3]";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_old_par0_1 $q_old_start $q_old_end $q_old_str > q_old_seq0_par_1`;
open(INPUT1,"q_old_seq0_par_1");
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
