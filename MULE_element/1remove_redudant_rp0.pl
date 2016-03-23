#!/usr/bin/perl -w
use strict;
my @file = ("MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new2.out.gff","MULE_gene/data/O.brachyantha.v1.4.fasta_new2.out.gff","MULE_gene/data/O.nivara_v1.0.pseudo.fasta_new2.out.gff","MULE_gene/data/O.glumaepatula.v1.5.fasta_new2.out.gff","MULE_gene/data/O.sativa.indica_9311.fasta_new2.out.gff","MULE_gene/data/Oryza_rufipogon_W1943_pseudo_v1.1.fasta_new2.out.gff","MULE_gene/data/O.barthii_v1.3.fasta_new2.out.gff","MULE_gene/data/Lperr_v1.4.fasta_new2.out.gff","MULE_gene/data/O.glaberrima.v1.fasta_new2.out.gff","MULE_gene/data/O.meridionalis.v1.3.fasta_new2.out.gff","MULE_gene/data/O.punctata_v1.2.fasta_new2.out.gff");
#my @file = ("MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new3.out.gff","MULE_gene/data/O.brachyantha.v1.4.fasta_new3.out.gff","MULE_gene/data/O.nivara_v1.0.pseudo.fasta_new3.out.gff","MULE_gene/data/O.glumaepatula.v1.5.fasta_new3.out.gff","MULE_gene/data/O.sativa.indica_9311.fasta_new3.out.gff","MULE_gene/data/Oryza_rufipogon_W1943_pseudo_v1.1.fasta_new3.out.gff","MULE_gene/data/O.barthii_v1.3.fasta_new3.out.gff","MULE_gene/data/Lperr_v1.4.fasta_new3.out.gff","MULE_gene/data/O.glaberrima.v1.fasta_new3.out.gff","MULE_gene/data/O.meridionalis.v1.3.fasta_new3.out.gff","MULE_gene/data/O.punctata_v1.2.fasta_new3.out.gff","MULE_gene/data/O.longistaminata.v0117.fasta_new3.out.gff");



for my $file (@file){
open(INPUT,"$file");
my $outfile = $file."_screen";
open(OUTPUT,">$outfile");
<INPUT>;
<INPUT>;
<INPUT>;
my $line = <INPUT>;
chomp $line;
my @info = split(/\s+/,$line);
my %hash;
my $prev_start = $info[3];
my $prev_end = $info[4];
my $prev_score = $info[5];
my $prev_str = $info[6];
my $prev_target = $info[9];
my $prev_target_start = $info[10];
my $prev_target_end = $info[11];
my $prev_chr = $info[0];
my $overlap=0;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $cur_start = $info[3];
my $cur_end = $info[4];
my $cur_score = $info[5];
my $cur_str = $info[6];
my $cur_target = $info[9];
my $cur_target_start = $info[10];
my $cur_target_end = $info[11];
my $cur_chr = $info[0];

	if(($prev_chr eq $cur_chr) && ($prev_str eq $cur_str)){
		my $temp_dist1 = abs($prev_start - $prev_end);
		my $temp_dist2 = abs($cur_start - $cur_end);
		my $overlap_len=0;
		if($temp_dist1 < $temp_dist2){
		$overlap_len = $temp_dist1/2;
		}else{
		$overlap_len = $temp_dist2/2;
		}
		$overlap=overlap($cur_start,$cur_end,$prev_start,$prev_end,$overlap_len);
		if($overlap ==1){
  		       if($prev_score > $cur_score){
			$cur_start = $prev_start;
			$cur_end = $prev_end;
			$cur_score = $prev_score;
			$cur_str = $prev_str;
			$cur_target = $prev_target;
			$cur_target_start =$prev_target_start;
			$cur_target_end = $prev_target_end;
			$cur_chr = $prev_chr;
			
		        }else{
			;
		        }
		}elsif($overlap ==0){
			if(exists $hash{$prev_chr}){
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{start} = $prev_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{end} = $prev_end;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{score} = $prev_score;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{str} = $prev_str;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target} = $prev_target;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_start} = $prev_target_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_end} = $prev_target_end;
				$hash{$prev_chr}->{count}++;
			}else{
			       $hash{$prev_chr}->{array}->[0]->{start} = $prev_start;
				$hash{$prev_chr}->{array}->[0]->{end} = $prev_end;
				$hash{$prev_chr}->{array}->[0]->{score} = $prev_score;
				$hash{$prev_chr}->{array}->[0]->{str} = $prev_str;
				$hash{$prev_chr}->{array}->[0]->{target} = $prev_target;
				$hash{$prev_chr}->{array}->[0]->{target_start} = $prev_target_start;
				$hash{$prev_chr}->{array}->[0]->{target_end} = $prev_target_end;
				$hash{$prev_chr}->{count}=1;
 			}
		}
	}else{
			if(exists $hash{$prev_chr}){
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{start} = $prev_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{end} = $prev_end;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{score} = $prev_score;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{str} = $prev_str;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target} = $prev_target;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_start} = $prev_target_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_end} = $prev_target_end;
				$hash{$prev_chr}->{count}++;
			}else{
			       $hash{$prev_chr}->{array}->[0]->{start} = $prev_start;
				$hash{$prev_chr}->{array}->[0]->{end} = $prev_end;
				$hash{$prev_chr}->{array}->[0]->{score} = $prev_score;
				$hash{$prev_chr}->{array}->[0]->{str} = $prev_str;
				$hash{$prev_chr}->{array}->[0]->{target} = $prev_target;
				$hash{$prev_chr}->{array}->[0]->{target_start} = $prev_target_start;
				$hash{$prev_chr}->{array}->[0]->{target_end} = $prev_target_end;
				$hash{$prev_chr}->{count}=1;
 			}
	}

	$prev_start = $cur_start;
	$prev_end = $cur_end;
	$prev_score = $cur_score;
	$prev_str = $cur_str;
	$prev_target = $cur_target;
	$prev_target_start = $cur_target_start;
	$prev_target_end = $cur_target_end;
	$prev_chr=$cur_chr;

}

 				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{start} = $prev_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{end} = $prev_end;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{score} = $prev_score;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{str} = $prev_str;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target} = $prev_target;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_start} = $prev_target_start;
				$hash{$prev_chr}->{array}->[$hash{$prev_chr}->{count}]->{target_end} = $prev_target_end;
				$hash{$prev_chr}->{count}++;
	

for my $key (sort {$a<=>$b} keys %hash){
		for(my $i=0; $i< $hash{$key}->{count};$i++){
			print OUTPUT "$key\t$hash{$key}->{array}->[$i]->{start}\t$hash{$key}->{array}->[$i]->{end}\t$hash{$key}->{array}->[$i]->{score}\t$hash{$key}->{array}->[$i]->{str}\t$hash{$key}->{array}->[$i]->{target}\t$hash{$key}->{array}->[$i]->{target_start}\t$hash{$key}->{array}->[$i]->{target_end}\n";
		}
}

}

sub overlap{
my ($start1,$end1,$start2,$end2,$overlap_len)= (@_);
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
goto line1;
}elsif($end1 >= $start2){
	if($start1 < $start2){
	  my $len = abs($end1 -$start2) +1;
	  if($len >=$overlap_len){
	  $overlap=1;
	  }
          goto line1;
        }elsif($start1 >= $start2){
            if($end1 < $end2){
	       my $len = abs($end1 - $start1)+1;
		if($len >= $overlap_len){
		$overlap=1;
		}
		goto line1;
	     }elsif($end1 >= $end2){
		if($start1 < $end2){
		  my $len = abs($end2 - $start1)+1;
		   if($len >=$overlap_len){
		   $overlap=1;
                   }
                   goto line1;
                }elsif($start1 >= $end2){
		   goto line1;
		}
	      }
        }
}

line1: return($overlap);
}






