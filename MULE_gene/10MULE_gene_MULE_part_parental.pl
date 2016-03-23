#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;

my @spe  =("japo");
for my $spe(@spe){
my $file = "MULE_gene/data/$spe"."_MULE_parental_seq_genome_blastn_gene_new1_exon30";
my $rank_file = "$file"."_rank";
my $output = $file."_parental";


open(OUTPUT,">$output");
my $mule_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5";
my $mule_gff1 = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new3.out.gff";
my $mule_gff = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new4.out.gff";

$rank_file = best_hit("$file","$rank_file");
open(INPUT,$rank_file);
my %count_mule=();
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $query_start = $info[6];
my $query_end = $info[7];
my $overlap_hit=0;
if(exists $count_mule{$info[0]}){
open(OUTPUT1,">hit_mule1");
for(my $i=0; $i<$count_mule{$info[0]}->{count}; $i++){
print OUTPUT1 "$count_mule{$info[0]}->{array}->[$i]->{start}\t$count_mule{$info[0]}->{array}->[$i]->{end}\n";

}
close(OUTPUT1);

$overlap_hit = overlap("chrg",$query_start,$query_end,"hit_mule1");
}
if($overlap_hit==0){
print "$line\n";
my $hit_chr = $info[1];
my $hit_start;
my $hit_end ;
my $hit_str;
if($info[8] <= $info[9] ){
$hit_str = "+";
$hit_start = $info[8];
$hit_end  = $info[9];
}elsif($info[9] < $info[8]){
$hit_str = "-";
$hit_start = $info[9];
$hit_end = $info[8];
}
`awk '{if(\$2 == "$hit_chr"){split(\$3,a,":");split(\$4,b,":"); c=a[1]-1000; d=b[2]+1000; print c,d;}}' $mule_file > chr_mule1`;

my $overlap_mule=overlap($hit_chr,$hit_start,$hit_end,"chr_mule1");

if($overlap_mule==1){
}elsif($overlap_mule==0){	
`awk '{if(\$1 =="$hit_chr"){a=\$4-1000; b=\$5+1000; print a,b;}}' $mule_gff > chr_mule1`;

my $overlap_mule=overlap($hit_chr,$hit_start,$hit_end,"chr_mule1");
if($overlap_mule==1){
;
}else{
  if(!(exists $count_mule{$info[0]}->{count})){
  $count_mule{$info[0]}->{count}=0;
  }
  $count_mule{$info[0]}->{array}->[$count_mule{$info[0]}->{count}]->{start} = $query_start;
  $count_mule{$info[0]}->{array}->[$count_mule{$info[0]}->{count}]->{end} = $query_end;
  $count_mule{$info[0]}->{count}++;
 
  $count_mule{$info[0]}->{hit} .= "$hit_chr:$hit_start:$hit_end:$hit_str\t";
}
}

}
}
for my $key (sort keys %count_mule){
print OUTPUT "$key\t$count_mule{$key}->{hit}\n";
}

}

sub overlap{

my ($cds_chr,$cds_start,$cds_end,$chr_file)=(@_);
open(INPUT1,"$chr_file");
my $exp_flag=0;
my $exp_str="";
my $overlap;
my $rna_info;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_start = $info1[0];
my $rna_end = $info1[1];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
	#	$overlap = abs($rna_end - $cds_start)/($cds_end - $cds_start+1);
		#if($overlap < 0.3){
	#	if($overlap <=0){

	#	$exp_flag=0;
	#	}else{
		$exp_flag=1;
		goto line2;

	#	}
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$overlap =1;
		goto line2;
		}elsif($rna_end <= $cds_start){
		$overlap=0;
		$exp_flag=0;
		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
	#	$overlap = abs($rna_end - $rna_start)/($cds_end - $cds_start +1);
#		if($overlap <0.3){
	#	if($overlap <=0){

	#	$exp_flag=0;
	#	}else{
		$exp_flag=1;
		goto line2;

	#	}

		}elsif(($rna_end >= $cds_end) &&($rna_start <= $cds_end)){
	#		$overlap = abs($cds_end - $rna_start)/($cds_end - $cds_start + 1);
#			if($overlap < 0.3){
	#		if($overlap <=0){

	#		$exp_flag=0;
	#		}else{
			$exp_flag=1;
			goto line2;

	#		}
		}elsif(($rna_end > $cds_end) && ($rna_start > $cds_end)){
			$exp_flag=0;
			
		}
	}
}

line2: return($exp_flag);
}


sub blat_score{
	my $psl = shift(@_);
	my @info = split(/\s+/,$psl);
	my $score = ($info[0] + $info[2]/2) - $info[1] - $info[4] - $info[6]; 
	return($score);
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
	$qAliSize = $sizeMul * ($info[12] - $info[11]+1);
	$tAliSize = $info[16] - $info[15]+1;
	my $qcover = ($info[0]+$info[2]/2)/$qAliSize;
	my $qcover1 = $qAliSize/$info[10];
	my $tcover = ($info[0] +$info[2]/2)/$tAliSize;

	if($qAliSize <= $tAliSize){
		$aliSize = $qAliSize;
	}else{
		$aliSize = $tAliSize;
	}

	if($aliSize <= 0){#this will not happen --ZCJ
		goto line1;
	}
	$sizeDif = $qAliSize - $tAliSize;
#	$sizeDif =0;
	if($sizeDif < 0){
		$sizeDif = 0;
	}
	$insertFactor = $info[4];
#	$insertFactor = 0;
#	$insertFactor += $info[6];
	$total = ($sizeMul * ($info[0] + $info[1]+$info[2]));
	if($total !=0){
		$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + round(3*log(1+$sizeDif))))/$total;
line1: return ($qcover1,$qcover,$tcover,$milliBad);
	}
}


sub round {
my $num = shift(@_);
if($num == (int($num)+0.5)){
return(int($num));
}else{
return(int($num+0.5));
}
}

sub best_hit{
my ($blat_out,$new_blat_out) = (@_);
open(OUTPUTb,">$new_blat_out");
open(INPUTb,$blat_out);
my %hash;
for(my $i=0; $i<=4; $i++){
my $temp_line = <INPUTb>;
chomp $temp_line;
}
my %count;
while(my $line = <INPUTb>){
if($line){
	chomp $line;
	my @info = split(/\s+/,$line);
	#my ($qcover1,$qcover,$tcover,$iden) = pslCal($line);
	my @name = split(/:/,$info[0]);
	my $qcover = $info[3]/($name[3] - $name[2]+1);
	my $qcover1 = $info[3] / ($info[7] - $info[6]+1);
        my $tcover;
	if($info[9] >= $info[8]){ 
	$tcover = $info[3] / ($info[9]- $info[8]+1);
	}else{
	$tcover = $info[3] / ($info[8]- $info[9]+1);
	}
	my $iden_score = $info[2];
	my $blat_score = $info[10];
	if(!(exists $count{$info[0]})){
	$count{$info[0]}=0;
	}
	$hash{$info[0]}->{$count{$info[0]}}->{iden_score} = $iden_score;
	$hash{$info[0]}->{$count{$info[0]}}->{score} = $blat_score;
	$hash{$info[0]}->{$count{$info[0]}}->{qcover} = $qcover;
	$hash{$info[0]}->{$count{$info[0]}}->{qcover1} = $qcover1;

	$hash{$info[0]}->{$count{$info[0]}}->{tcover} = $tcover;
	$hash{$info[0]}->{$count{$info[0]}}->{line} = $line;
	$count{$info[0]}++;
}
}
my %count_key;
for my $key (sort keys %hash){
	my %sub_hash = %{$hash{$key}};
	for my $sub_key (sort {	$sub_hash{$a}->{score} <=> $sub_hash{$b}->{score}
	||
	$sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score} } keys %sub_hash){
		if($sub_hash{$sub_key}->{score} < 1e-10){
		print OUTPUTb "$sub_hash{$sub_key}->{line}\t$sub_hash{$sub_key}->{qcover1}\t$sub_hash{$sub_key}->{qcover}\t$sub_hash{$sub_key}->{tcover}\n";
		}
	}
}
close OUTPUTb;
return("$new_blat_out");
}


