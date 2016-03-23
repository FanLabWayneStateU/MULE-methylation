#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;

my $spe = $ARGV[0];
my $blat_spe = $ARGV[1];
my $genome_seq_query = "$ARGV[2]"; #"MULE_gene/data/O.brachyantha.v1.4.fasta";
my $genome_seq_target = "$ARGV[3]"; #"MULE_gene/data/O.nivara_v1.0.pseudo.fasta";

my %query_hash = %{get_seq($genome_seq_query)};
my %target_hash = %{get_seq($genome_seq_target)};


#for my $spe (@spe){
#my $blat_out = "MULE_gene/data/MULE_$spe"."_flank_genome_niva.psl";
my $blat_out = "$ARGV[4]";
my $blat_out_best = $blat_out."_best";
$blat_out_best=best_hit("$blat_out","$blat_out_best");

#my $flank_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_no_TE60or_flank_list";
my $flank_file = "$ARGV[5]";
#my $output = $flank_file."_niva_ortho1";
my $output = $flank_file."_$blat_spe"."_ortho1";
open(OUTPUT,">$output");

open(INPUT,"$flank_file");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $query = $info[0];
my ($query_a1,$query_a2,$query_b1,$query_b2) = ($info[2],$info[1],$info[3],$info[4]);
my $query_ortho = 0; 
my $query_ortho_best=0;
$query_ortho = judge_ortho($query,$query_a1,$query_a2,$query_b1,$query_b2,$blat_out_best,\%query_hash,\%target_hash,$spe);
my ($query_ortho_cover,$query_blat_chr,$query_blat_start,$query_blat_end,$query_blat_str) = cover($query,$blat_out_best);
my ($query_ortho_cover1,$query_blat_chr1,$query_blat_start1,$query_blat_end1,$query_blat_str1) = cover($query_a1,$blat_out_best);
my ($query_ortho_cover2,$query_blat_chr2,$query_blat_start2,$query_blat_end2,$query_blat_str2) = cover($query_a2,$blat_out_best);
my ($query_ortho_cover3,$query_blat_chr3,$query_blat_start3,$query_blat_end3,$query_blat_str3) = cover($query_b1,$blat_out_best);
my ($query_ortho_cover4,$query_blat_chr4,$query_blat_start4,$query_blat_end4,$query_blat_str4) = cover($query_b2,$blat_out_best);

if(($query_ortho_cover >= 0.30) && ($query_ortho==1)){
$query_ortho_best=1;
}
print OUTPUT "$query\t$query_ortho_best\t$query_ortho_cover:$query_blat_chr:$query_blat_start:$query_blat_end:$query_blat_str\t$query_ortho_cover1:$query_blat_chr1:$query_blat_start1:$query_blat_end1:$query_blat_str1\t$query_ortho_cover2:$query_blat_chr2:$query_blat_start2:$query_blat_end2:$query_blat_str2\t$query_ortho_cover3:$query_blat_chr3:$query_blat_start3:$query_blat_end3:$query_blat_str3\t$query_ortho_cover4:$query_blat_chr4:$query_blat_start4:$query_blat_end4:$query_blat_str4\n";
}
#}
sub cover{
my ($gene_id,$blat_best) = (@_);
`grep "$gene_id" $blat_best | cat > gene_ortho_cov101_perr2`;
open(INPUTc,"gene_ortho_cov101_perr2");
my $linec = <INPUTc>;
if($linec){
chomp $linec;
my @infoc = split(/\s+/,$linec);
my $cov=0;
if($infoc[10] >0){
$cov = (abs($infoc[12] - $infoc[11])) / $infoc[10];
}

#print "$cov,$infoc[13],$infoc[15],$infoc[16],$infoc[8]\n";
#exit;
return($cov,$infoc[13],$infoc[15],$infoc[16],$infoc[8]);
}else{
return(0,0,0,0,0);
}
}

sub judge_ortho{
my ($query,$query_a1,$query_a2,$query_b1,$query_b2,$blat_out_best,$query_hash,$target_hash,$spe) = (@_);
my %query_hash = %{$query_hash};
my %target_hash = %{$target_hash};

my $ortho_count=0;
my ($q_ortho_chr,$q_ortho_start,$q_ortho_end,$q_ortho_str) = get_ortho($query,$blat_out_best);
my ($q_ortho_chr_a1,$q_ortho_start_a1,$q_ortho_end_a1,$q_ortho_str_a1) = get_ortho($query_a1,$blat_out_best);
my ($q_ortho_chr_a2,$q_ortho_start_a2,$q_ortho_end_a2,$q_ortho_str_a2) = get_ortho($query_a2,$blat_out_best);
my ($q_ortho_chr_b1,$q_ortho_start_b1,$q_ortho_end_b1,$q_ortho_str_b1) = get_ortho($query_b1,$blat_out_best);
my ($q_ortho_chr_b2,$q_ortho_start_b2,$q_ortho_end_b2,$q_ortho_str_b2) = get_ortho($query_b2,$blat_out_best);
my $new_q_ortho_start = $q_ortho_start - 4000;
my $new_q_ortho_end = $q_ortho_end + 4000;
my @name = split(/:/,$query);
my $query_chr;
if($name[1]=~ /(\d+)/){
$query_chr = $1;
if($query_chr =~ /^0(\d+)/){
$query_chr = $1;
}
}
my $query_len = $name[3]-$name[2]+1;
#my @name = split(/:/,$query);
my $query_chr_ortho;
if($q_ortho_chr=~ /(\d+)/){
$query_chr_ortho = $1;
if($query_chr_ortho =~ /^0(\d+)/){
$query_chr_ortho = $1;
}
}

if($query_chr_ortho == $query_chr){
$ortho_count++;
}
#if(($q_ortho_chr_a1 eq $q_ortho_chr) && (($q_ortho_start - 4000) <= $q_ortho_start_a1) && (($q_ortho_end + 4000) >= $q_ortho_end_a1)){
if($q_ortho_chr_a1 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_a1,$q_ortho_end_a1,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){
$ortho_count++;
}
}
#if(($q_ortho_chr_a2 eq $q_ortho_chr) && (($q_ortho_start - 4000) <= $q_ortho_start_a2) && (($q_ortho_end + 4000) >= $q_ortho_end_a2)){
if($q_ortho_chr_a2 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_a2,$q_ortho_end_a2,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){
$ortho_count++;
}
}
#if(($q_ortho_chr_b1 eq $q_ortho_chr) && (($q_ortho_start - 4000) <= $q_ortho_start_b1) && (($q_ortho_end + 4000) >= $q_ortho_end_b1)){
if($q_ortho_chr_b1 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_b1,$q_ortho_end_b1,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){

$ortho_count++;
}
}
#if(($q_ortho_chr_b2 eq $q_ortho_chr) && (($q_ortho_start - 4000) <= $q_ortho_start_b2) && (($q_ortho_end + 4000) >= $q_ortho_end_b2)){
if($q_ortho_chr_b2 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_b2,$q_ortho_end_b2,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){


$ortho_count++;
}
}
if($query_chr_ortho != $query_chr){
my $q_ortho_chr_a2_new = get_chr_id($q_ortho_chr_a2);
my $q_ortho_chr_b1_new = get_chr_id($q_ortho_chr_b1);
if(($q_ortho_chr_a2_new == $q_ortho_chr_b1_new) && ($q_ortho_chr_a2_new == $query_chr) &&($q_ortho_str_a2 eq $q_ortho_str_b1)){
my $ortho_flank_diff = -1;
if(($q_ortho_start_a2 <= $q_ortho_start_b1) &&($q_ortho_end_a2 <= $q_ortho_end_b1)){
$ortho_flank_diff = $q_ortho_start_b1 - $q_ortho_end_a2 + 1;
}elsif(($q_ortho_start_a2 >= $q_ortho_start_b1) &&($q_ortho_end_a2 >= $q_ortho_end_b1)){
$ortho_flank_diff = $q_ortho_start_a2 - $q_ortho_end_b1 + 1;
}
if($ortho_flank_diff > 0 ){
my $ortho_insert_diff = abs($query_len - $ortho_flank_diff)+1;
my $ortho_insert_ratio = $ortho_insert_diff/$ortho_flank_diff;
if($query_len > $ortho_flank_diff){
$ortho_insert_ratio = $ortho_insert_diff/$query_len;
}
if($ortho_insert_ratio <= 0.7){
my ($left_tir_start,$left_tir_end,$right_tir_start,$right_tir_end)=get_tir($query,$spe);
my $left_tir_len = $left_tir_end -$left_tir_start +1;
my $left_tir_ortho_len = $q_ortho_end_a2 - $q_ortho_start_a2 + 1;
my $right_tir_len = $right_tir_end -$right_tir_start +1;
my $right_tir_ortho_len = $q_ortho_end_b1 - $q_ortho_start_b1 + 1;
my ($new_left_tir_start,$new_left_tir_end,$new_right_tir_start,$new_right_tir_end,$new_q_ortho_start_a2,$new_q_ortho_start_b1,$new_q_ortho_end_a2,$new_q_ortho_end_b1)= ($left_tir_start,$left_tir_end,$right_tir_start,$right_tir_end,$q_ortho_start_a2,$q_ortho_start_b1,$q_ortho_end_a2,$q_ortho_end_b1);

if(($q_ortho_start_a2 <= $q_ortho_start_b1) &&($q_ortho_end_a2 <= $q_ortho_end_b1)){

if($left_tir_len >= $left_tir_ortho_len){
$new_left_tir_start = ($left_tir_end - $left_tir_ortho_len);
}elsif($left_tir_len < $left_tir_ortho_len){
$new_q_ortho_start_a2 = ($new_q_ortho_end_a2 - $left_tir_len);
}

if($right_tir_len >= $right_tir_ortho_len){
$new_right_tir_end = ($right_tir_start + $right_tir_ortho_len);
}elsif($right_tir_len < $right_tir_ortho_len){
$new_q_ortho_end_b1 = ($new_q_ortho_start_b1 + $right_tir_len);
}

}elsif(($q_ortho_start_a2 >= $q_ortho_start_b1) &&($q_ortho_end_a2 >= $q_ortho_end_b1)){
if($left_tir_len >= $left_tir_ortho_len){
$new_left_tir_start = ($left_tir_end - $left_tir_ortho_len);
}elsif($left_tir_len < $left_tir_ortho_len){
$new_q_ortho_end_a2 = ($new_q_ortho_start_a2 + $left_tir_len);
}
if($right_tir_len >= $right_tir_ortho_len){
$new_right_tir_end = ($right_tir_start + $right_tir_ortho_len);
}elsif($right_tir_len < $right_tir_ortho_len){
$new_q_ortho_start_b1 = ($new_q_ortho_end_b1 - $right_tir_len);
}
}
my $left_tir_seq = get_seq_seg(\%query_hash,$new_left_tir_start,$new_left_tir_end,"+",$name[1]);
open(OUTPUTl,">left_tir_perr2");
print OUTPUTl ">left_tir_perr2\n";
print OUTPUTl "$left_tir_seq\n";

my $left_tir_seq_ortho = get_seq_seg(\%target_hash,$new_q_ortho_start_a2,$new_q_ortho_end_a2,$q_ortho_str_a2,$q_ortho_chr_a2);
open(OUTPUTl1,">left_tir_ortho_perr2");
print OUTPUTl1 ">left_tir_ortho_perr2\n";
print OUTPUTl1 "$left_tir_seq_ortho\n";
`bin/x86_64-redhat-linux-gnu/blat left_tir_perr2 left_tir_ortho_perr2 -minIdentity=30 left_tir_perr2.psl`;
my ($tir_cover_l,$tir_chr_l,$tir_start_l,$tir_end_l,$tir_str_l) = cover("left_tir_ortho_perr2","left_tir_perr2.psl");

my $right_tir_seq = get_seq_seg(\%query_hash,$new_right_tir_start,$new_right_tir_end,"+",$name[1]);
open(OUTPUTr,">right_tir_perr2");
print OUTPUTr ">right_tir_perr2\n";
print OUTPUTr "$right_tir_seq\n";

my $right_tir_seq_ortho = get_seq_seg(\%target_hash,$new_q_ortho_start_b1,$new_q_ortho_end_b1,$q_ortho_str_b1,$q_ortho_chr_b1);
open(OUTPUTr1,">right_tir_ortho_perr2");
print OUTPUTr1 ">right_tir_ortho_perr2\n";
print OUTPUTr1 "$right_tir_seq_ortho\n";
`bin/x86_64-redhat-linux-gnu/blat right_tir_perr2 right_tir_ortho_perr2 -minIdentity=30 right_tir_perr2.psl`;
#exit;
my ($tir_cover_r,$tir_chr_r,$tir_start_r,$tir_end_r,$tir_str_r) = cover("right_tir_ortho_perr2","right_tir_perr2.psl");
if(($tir_cover_l >= 0.5 ) &&($tir_cover_r >=0.5)){
$ortho_count=2;
}
}
}
}
}
my $ortho=0;
if($ortho_count >= 2){
$ortho=1;
}
return($ortho);
}

sub get_seq_seg{
my ($query_hash,$new_left_tir_start,$new_left_tir_end,$str,$chr)=(@_);
my %query_hash = %{$query_hash};
open(OUTPUT1,">scaffold_new_perr2");
print OUTPUT1 ">$chr\n";
print OUTPUT1 "$query_hash{$chr}\n";

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold_new_perr2 $new_left_tir_start $new_left_tir_end $str > middle_perr2`;
open(INPUTlf,"middle_perr2");
my $seq = <INPUTlf>;
chomp $seq;
return($seq);
}

sub get_ortho{
my ($gene_id,$blat_best) = (@_);
#print "$gene_id,$blat_best\n";
`grep "$gene_id" $blat_best | cat > gene_ortho_best101_perr2`;
open(INPUTb1,"gene_ortho_best101_perr2");
my $lineb = <INPUTb1>;
if($lineb){
chomp $lineb;
my @infob = split(/\s+/,$lineb);
return($infob[13],$infob[15],$infob[16],$infob[8]);
}else{
return(0,0,0)
}
}

sub get_tir{
my ($name,$spe) =(@_);
my @name = split(/:/,$name);
#my $tir_file = "MULE_gene/data/pack_MULE_candidate_"."$spe"."_jiang_repeatscout_no_TE60or";
my $tir_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

`grep "$name[0]" $tir_file | grep "$name[1]" | grep "$name[2]" | grep "$name[3]" > coor_line_perr2`;
open(INPUTc,"coor_line_perr2");
my $linec = <INPUTc>;
chomp $linec;
my @info = split(/\s+/,$linec);
my @left = split(/:/,$info[2]);
my @right = split(/:/,$info[3]);
return($left[0],$left[1],$right[0],$right[1]);
}
sub get_chr_id{
my $q_ortho_chr = shift(@_);
my $query_chr_ortho;
if($q_ortho_chr=~ /(\d+)/){
$query_chr_ortho = $1;
if($query_chr_ortho =~ /^0(\d+)/){
$query_chr_ortho = $1;
}
}
return($query_chr_ortho);
}

sub overlap{
my ($start1,$end1,$start2,$end2)= (@_);
if($start1 > $end1){
my $temp =$start1;
$start1 = $end1;
$end1 = $temp;
}
if($start2 > $end2){
my $temp = $start2;
$start2 = $end2;
$end2 = $temp;
}
my $overlap=0;
if($end1 < $start2){
goto line2;
}elsif($end1 >= $start2){
	if($start1 < $start2){
	  if($end1 > $end2){
	  my $len = abs($end2 - $start2)+1;
	  if($len >=20){
	  $overlap=1;
	  }
	  goto line2;
	  }elsif($end1 <=$end2){
	  my $len = abs($end1 -$start2) +1;
	  if($len >=20){
	  $overlap=1;
	  }
	  goto line2;
	  }
        }elsif($start1 >= $start2){
            if($end1 < $end2){
	       my $len = abs($end1 - $start1)+1;
		if($len >= 20){
		$overlap=1;
		}
		goto line2;
	     }elsif($end1 >= $end2){
		if($start1 < $end2){
		  my $len = abs($end2 - $start1)+1;
		   if($len >=20){
		   $overlap=1;
                   }
                   goto line2;
                }elsif($start1 >= $end2){
		   goto line2;
		}
	      }
        }
}

line2: return($overlap);
}



sub best_hit{
my ($blat_out,$new_blat_out) = (@_);
open(OUTPUTb,">$new_blat_out");
open(INPUTb,$blat_out);
my %hash;
for(my $i=0; $i<=4; $i++){
my $temp_line = <INPUTb>;
chomp $temp_line;
#print  OUTPUTb "$temp_line\n";
}
my %count;
while(my $line = <INPUTb>){
if($line){
#	print $line;
	#exit;
	chomp $line;
	my @info = split(/\s+/,$line);
#	print "line:$line\n";
	my ($qcover,$tcover,$iden) = pslCal($line);
	my $iden_score = 100.0 - $iden * 0.1;
	my $blat_score = blat_score($line);
	if(!(exists $count{$info[9]})){
	$count{$info[9]}=0;
	}
	$hash{$info[9]}->{$count{$info[9]}}->{iden_score} = $iden_score;
	$hash{$info[9]}->{$count{$info[9]}}->{score} = $blat_score;
	$hash{$info[9]}->{$count{$info[9]}}->{qcover} = $qcover;
	$hash{$info[9]}->{$count{$info[9]}}->{tcover} = $tcover;
	$hash{$info[9]}->{$count{$info[9]}}->{line} = $line;
	$count{$info[9]}++;
}
}
my %count_key;
for my $key (sort keys %hash){
	my %sub_hash = %{$hash{$key}};
#	sub by_score {
#	$sub_hash{$b}->{score} <=> $sub_hash{$a}->{score}
#	||
#	$sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score}
#	}
	for my $sub_key (sort {	$sub_hash{$b}->{score} <=> $sub_hash{$a}->{score}
	||
	$sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score} } keys %sub_hash){
		if((!(exists $count_key{$key}))){
		print OUTPUTb "$sub_hash{$sub_key}->{line}\n";
#		exit;
		$count_key{$key} =1;
		}
	}
}
close OUTPUTb;
#open(INPUTtest,"new_blatout");
#while(my $linet = <INPUTtest>){
#print "$linet";
#}
return("$new_blat_out");
}



sub blat_score{
	my $psl = shift(@_);
	my @info = split(/\s+/,$psl);
#	my $score = ($info[0] + $info[2]/2) - $info[1]; # - $info[5] - $info[7];
	my $score = ($info[0] + $info[2]/2) - $info[1]  - $info[4] - $info[6];

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
		$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + round(3*log(1+$sizeDif))))/$total;##where dose this come from? shall we cite paper?
line1: return ($qcover,$tcover,$milliBad);
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
