#!/usr/bin/perl -w
use strict;
my @spe = ("japo");
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";
my $output_total1 = $file."_tir_mafft_result";
open(OUTPUTtotal1,">$output_total1");
my $output_total2 = $file."_tir";
open(OUTPUTtotal2,">$output_total2");
my $genome_seq = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq_hash = %{get_seq($genome_seq)};

open(INPUT,$file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @left = split(/:/,$info[2]);
my @right = split(/:/,$info[3]);
my $chr=$info[1];
my ($left_start,$left_end,$left_str)=($left[0],$left[1],$left[2]);
my ($right_start,$right_end,$right_str)=($right[0],$right[1],$right[2]);
my $mule_name = "$info[0]:$info[1]:$left[1]:$right[0]";
my $left_name = "$mule_name:L";
my $right_name  = "$mule_name:R";
open(OUTPUT,">scaffold_new_new");
print OUTPUT ">$chr\n$seq_hash{$chr}\n";
open(OUTPUT1,">left_right");
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new_new $left_start $left_end $left_str > new_seq`;
open(INPUT2,"new_seq");
my $gene_seq1 = <INPUT2>;
chomp $gene_seq1;
print OUTPUT1 ">$left_name\n";
print OUTPUT1 "$gene_seq1\n";
print OUTPUTtotal2 ">$left_name\n";
print OUTPUTtotal2 "$gene_seq1\n";

`rm new_seq`;
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new_new $right_start $right_end $right_str > new_seq`;
open(INPUT3,"new_seq");
my $gene_seq2 = <INPUT3>;
chomp $gene_seq2;
print OUTPUT1 ">$right_name\n";
print OUTPUT1 "$gene_seq2\n";
print OUTPUTtotal2 ">$right_name\n";
print OUTPUTtotal2 "$gene_seq2\n";


`rm new_seq`;
`mafft-linux64/mafft.bat --localpair --maxiterate 100 left_right > left_right.out`;
my %align_seq = %{get_seq("left_right.out")};
my $left_seq = $align_seq{$left_name};
my $right_seq = $align_seq{$right_name};
my @left_seq = split(//,$left_seq);
my @right_seq = split(//,$right_seq);
my $left_length = $left_end - $left_start + 1;
my $right_length = $right_end - $right_start + 1;

my ($left_match,$right_match,$left_mis,$right_mis,$left_gap,$right_gap);
for(my $i=0; $i<=$#left_seq; $i++){
if($left_seq[$i] eq "-"){
$left_gap++;
}
if($right_seq[$i] eq "-"){
$right_gap++;
}
if(($left_seq[$i] ne "-")&&($right_seq[$i] ne "-")&& ($left_seq[$i] ne $right_seq[$i])){
$left_mis++;
$right_mis++;
}

if($left_seq[$i] eq $right_seq[$i]){
$left_match++;
$right_match++;
}

}


my $left_match1 = $left_match/$left_length;
my $right_match1 = $right_match/$right_length;
my $left_mis1 = $left_mis/$left_length;
my $right_mis1 = $right_mis/$right_length;
my $left_gap1 = $left_gap/$left_length;
my $right_gap1 = $right_gap/$right_length;
my $identity= $left_match/length($left_seq);

print OUTPUTtotal1 "$mule_name\t$left_length:$left_match:$left_mis:$left_gap\t$left_match1:$left_mis1:$left_gap1\t$right_length:$right_match:$right_mis:$right_gap\t$right_match1:$right_mis1:$right_gap1\t$identity\n";

}
}



sub pslCal{
	my $psl = shift(@_);
	my $sizeMul = 1;
	my ($qAliSize,$tAliSize,$aliSize);
	my $milliBad=0;
	my $sizeDif;
	my $insertFactor; 
	my $total;
	my @info = split(/\s+/,$psl);
	$qAliSize = $sizeMul * ($info[12] - $info[11]);
	$tAliSize = $info[16] - $info[15];
	my $qcover = $qAliSize/$info[10];
	my $tcover = $tAliSize/$info[14];
	if($qAliSize <= $tAliSize){
		$aliSize = $qAliSize;
	}else{
		$aliSize = $tAliSize;
	}

	if($aliSize <= 0){#this will not happen --ZCJ
		goto line1;
	}
	$sizeDif = $qAliSize - $tAliSize;
	if($sizeDif < 0){
#		$sizeDif = -$sizeDif;
		$sizeDif =0;
	}
	
	#$insertFactor = $info[5];
	$insertFactor = $info[4];

#	$insertFactor += $info[7];
	$total = ($sizeMul * ($info[0] + $info[1]+$info[2]));
	if($total !=0){
		$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + ceil(3*log(1+$sizeDif))))/$total;
line1: return ($qcover,$tcover,$milliBad);
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
