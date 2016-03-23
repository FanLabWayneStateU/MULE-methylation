#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;

my $blat_out = "$ARGV[0]";
my $blat_out_best = $blat_out."_best_sort";
$blat_out_best=best_hit("$blat_out","$blat_out_best");

sub cover{
my ($gene_id,$blat_best,$temp_blat) = (@_);
`awk '{if(\$10=="$gene_id"){print}}' $blat_best | cat > $temp_blat`;
open(INPUTc,"$temp_blat");
my $linec = <INPUTc>;
if($linec){
chomp $linec;
my @infoc = split(/\s+/,$linec);
my $cov=0;
if($infoc[10] >0){
$cov = (abs($infoc[12] - $infoc[11])) / $infoc[10];
}

return($cov,$infoc[13],$infoc[15],$infoc[16],$infoc[8]);
}else{
return(0,0,0,0,0);
}
}

sub get_ortho{
my ($gene_id,$blat_best,$temp_blat) = (@_);
`awk '{if(\$10=="$gene_id"){print}}' $blat_best | cat > $temp_blat`;

open(INPUTb1,"$temp_blat");
my $lineb = <INPUTb1>;
if($lineb){
chomp $lineb;
my @infob = split(/\s+/,$lineb);
return($infob[13],$infob[15],$infob[16],$infob[8]);
}else{
return(0,0,0,0)
}
}

sub judge_ortho{
my ($query,$query_a1,$query_a2,$query_b1,$query_b2,$blat_out_best,$temp_blat) = (@_);

my $ortho_count=0;
my ($q_ortho_chr,$q_ortho_start,$q_ortho_end,$q_ortho_str) = get_ortho($query,$blat_out_best,$temp_blat);
my ($q_ortho_chr_a1,$q_ortho_start_a1,$q_ortho_end_a1,$q_ortho_str_a1) = get_ortho($query_a1,$blat_out_best,$temp_blat);
my ($q_ortho_chr_a2,$q_ortho_start_a2,$q_ortho_end_a2,$q_ortho_str_a2) = get_ortho($query_a2,$blat_out_best,$temp_blat);
my ($q_ortho_chr_b1,$q_ortho_start_b1,$q_ortho_end_b1,$q_ortho_str_b1) = get_ortho($query_b1,$blat_out_best,$temp_blat);
my ($q_ortho_chr_b2,$q_ortho_start_b2,$q_ortho_end_b2,$q_ortho_str_b2) = get_ortho($query_b2,$blat_out_best,$temp_blat);
my $new_q_ortho_start = $q_ortho_start - 4000;
my $new_q_ortho_end = $q_ortho_end + 4000;
my $query_chr;
my $query_chr_ortho;
if($q_ortho_chr_a1 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_a1,$q_ortho_end_a1,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){
$ortho_count++;
}
}
if($q_ortho_chr_a2 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_a2,$q_ortho_end_a2,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){
$ortho_count++;
}
}
if($q_ortho_chr_b1 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_b1,$q_ortho_end_b1,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){

$ortho_count++;
}
}
if($q_ortho_chr_b2 eq $q_ortho_chr){
my $overlap = overlap($q_ortho_start_b2,$q_ortho_end_b2,$new_q_ortho_start,$new_q_ortho_end);
if($overlap==1){


$ortho_count++;
}
}

my $ortho=0;
if($ortho_count >= 2){
$ortho=1;
}
return($ortho);
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
	  if($len >=50){
	  $overlap=1;
	  }
	  goto line2;
	  }elsif($end1 <=$end2){
	  my $len = abs($end1 -$start2) +1;
	  if($len >=50){
	  $overlap=1;
	  }
	  goto line2;
	  }
        }elsif($start1 >= $start2){
            if($end1 < $end2){
	       my $len = abs($end1 - $start1)+1;
		if($len >= 50){
		$overlap=1;
		}
		goto line2;
	     }elsif($end1 >= $end2){
		if($start1 < $end2){
		  my $len = abs($end2 - $start1)+1;
		   if($len >=50){
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
}
my %count;
while(my $line = <INPUTb>){
if($line){
	chomp $line;
	my @info = split(/\s+/,$line);
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
	$hash{$info[9]}->{$count{$info[9]}}->{target} = $info[13];
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
		my @query_id = split(/:/,$key);
		my $query_id = "$query_id[0]:$query_id[1]:$query_id[2]:$query_id[3]:$query_id[4]";
		my @target_id = split(/:/,$sub_hash{$sub_key}->{target});
		my $target_id = "$target_id[0]:$target_id[1]:$target_id[2]:$target_id[3]:$target_id[4]";
		#if(($sub_hash{$sub_key}->{iden_score} >=50)){

#		if(($sub_hash{$sub_key}->{iden_score} >=50)&&(($sub_hash{$sub_key}->{qcover}>=.7)||($sub_hash{$sub_key}->{tcover}>=.7))){
		if(($sub_hash{$sub_key}->{iden_score} >=50)&&(($sub_hash{$sub_key}->{qcover}>=.30)||($sub_hash{$sub_key}->{tcover}>=.30))){
		if(($query_id ne $target_id)&&($query_id[0] eq $target_id[0])){
		my @info = split(/\s+/,$sub_hash{$sub_key}->{line});
		print OUTPUTb "$info[9]\t$info[10]\t$info[11]\t$info[12]\t$info[13]\t$info[14]\t$info[15]\t$info[16]\t$sub_hash{$sub_key}->{iden_score}\t$sub_hash{$sub_key}->{qcover}\t$sub_hash{$sub_key}->{tcover}\t$info[8]\n";
		
#		exit;
		$count_key{$key} =1;
		}
		}
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
	#print "$psl";
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
		$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + round(3*log(1+$sizeDif))))/$total;
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
