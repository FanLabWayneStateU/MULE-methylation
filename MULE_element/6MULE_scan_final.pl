#!/usr/bin/perl -w
use strict;

my $TE_file = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new3.out.gff_screen";

my $output = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5";

my $genome_seq = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";


my %seq_hash = %{get_seq($genome_seq)};
open(OUTPUT1,">$output");
my %hash;

open(INPUT,"$TE_file");
while(my $line = <INPUT>){
chomp $line;
my @info1 = split(/\s+/,$line);
if(!exists $hash{$info1[0]}){
$hash{$info1[0]}->{count}=0;
}
my @te_name = split(/[:|\"]/,$info1[5]);
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{start} = $info1[1];
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{end} = $info1[2];
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{str} = $info1[4];
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{te_start} = $info1[6];
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{te_end} = $info1[7];
$hash{$info1[0]}->{array}->[$hash{$info1[0]}->{count}]->{name} = $te_name[2];

$hash{$info1[0]}->{count}++;
}
for my $key (sort {$a cmp $b} keys %hash){

	my %mark;
	open(OUTPUT,">scaffold2_new2_2");
	print OUTPUT ">$key\n";
	print OUTPUT "$seq_hash{$key}\n";


for (my $i=0; $i<($hash{$key}->{count}-1); $i++){
		if(exists $mark{$i}){
		;
		}else{
				for (my $j=$i+1; $j <= $i+1;$j++){


			if(exists $mark{$j}){
			;
			}else{

			if(($hash{$key}->{array}->[$i]->{name} =~ /(\S+)L$/)){
			$hash{$key}->{array}->[$i]->{name}="$1";
			}
  			if(($hash{$key}->{array}->[$j]->{name} =~ /(\S+)R$/)){

			$hash{$key}->{array}->[$j]->{name}="$1";
			}
			if((($hash{$key}->{array}->[$i]->{name} =~ /L$/) && ($hash{$key}->{array}->[$i]->{name} =~ /L$/))||(($hash{$key}->{array}->[$i]->{name} =~ /R$/) && ($hash{$key}->{array}->[$i]->{name} =~ /R$/))){
			$mark{$i}=1;
			$mark{$j}=1;
			goto line1;
			}


	 		if(($hash{$key}->{array}->[$i]->{str} eq "+") && ( $hash{$key}->{array}->[$j]->{str} eq "-") && ($hash{$key}->{array}->[$i]->{name} eq $hash{$key}->{array}->[$j]->{name})){
			my $dist = $hash{$key}->{array}->[$j]->{start} - $hash{$key}->{array}->[$i]->{end} +1;
		#	if($dist <=20000){			
			if($dist <=40000){				
	
		#			my $overlap=1;

		#			if( $overlap==1){
	
		my $origin_start = $hash{$key}->{array}->[$i]->{start};
		my $TSD_start1 = $origin_start-21;
		my $TSD_end1 = $origin_start-1;
		`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold2_new2_2 $TSD_start1 $TSD_end1 + > 5_TSD_2_new2_2`;
					open(INPUT3,"5_TSD_2_new2_2");
					my $line3 = <INPUT3>;
					chomp $line3;
					$line3 =uc($line3);
					open(OUTPUT3,">5_TSD2_1_new2_2");
					my $num=0; 
					my %hash_seq;
					for(my $temp_start = 0; $temp_start <= 10; $temp_start++){
					#my $temp_start = 10;
					my @len = (11,10,9,8,7);
					#my @len = (11);
					for my $len (@len){
                                        my $temp_seq = substr($line3,$temp_start,$len);
                                        my $temp_name = "$temp_start"."_$len";
					print OUTPUT3 ">$temp_name\n";
					print OUTPUT3 "$temp_seq\n";
					$hash_seq{$temp_name} = $temp_seq;
					}
					}
					my @len = (10,9,8,7);
					#my @len = (10);
					my $temp_start =11;
					for my $len (@len){
					my $temp_seq = substr($line3,11,$len);
                                        my $temp_name = "$temp_start"."_$len";
					print OUTPUT3 ">$temp_name\n";
					print OUTPUT3 "$temp_seq\n";
					$hash_seq{$temp_name} = $temp_seq;
					}
					my @len = (9,8,7);
					#my @len = (9);
					my $temp_start = 12;
					for my $len (@len){
					my $temp_seq = substr($line3,12,$len);
                                        my $temp_name = "$temp_start"."_$len";
					print OUTPUT3 ">$temp_name\n";
					print OUTPUT3 "$temp_seq\n";
					$hash_seq{$temp_name} = $temp_seq;
					}
					my @len = (8,7);
					#my @len = (8);
					my $temp_start = 13;
					for my $len (@len){
					my $temp_seq = substr($line3,13,$len);
                                        my $temp_name = "$temp_start"."_$len";
					print OUTPUT3 ">$temp_name\n";
					print OUTPUT3 "$temp_seq\n";
					$hash_seq{$temp_name} = $temp_seq;
					}
					my @len = (7);
					my $temp_start =14;
					for my $len (@len){
					my $temp_seq = substr($line3,14,$len);
                                        my $temp_name = "$temp_start"."_$len";
					print OUTPUT3 ">$temp_name\n";
					print OUTPUT3 "$temp_seq\n";
					$hash_seq{$temp_name} = $temp_seq;
					}



					my $origin_end = $hash{$key}->{array}->[$j]->{end};

					my $TSD_start2 = $origin_end+1;
					my $TSD_end2 = $origin_end+21;
					#my $TSD_end2 = $origin_end+11;

					`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl scaffold2_new2_2 $TSD_start2 $TSD_end2 + > 3_TSD_2_new2_2`;
					open(INPUT4,"3_TSD_2_new2_2");
					my $line4 = <INPUT4>;
					chomp $line4;
					$line4 = uc($line4);
					open(OUTPUT4,">3_TSD2_1_new2_2");
					print OUTPUT4 ">3\n";
					print OUTPUT4 "$line4";
					`bin/patman -D 3_TSD2_1_new2_2 -P 5_TSD2_1_new2_2 -s -g 2 -o TSD2_new2_2.test`;

					 open(INPUT2,"TSD2_new2_2.test");
					 my %result;
					 while(my $line2 = <INPUT2>){
					 chomp $line2;
					 my @info2 = split(/\s+/,$line2);
                                         my @pattern = split(/\_/,$info2[1]);
					 my $id = "$info2[1]"."_$info2[2]"."_$info2[5]";
					 $result{$id}->{id} = $info2[1];	

                                         $result{$id}->{pattern_len} = $pattern[1];
                                         $result{$id}->{mis} = $info2[5];
                                         $result{$id}->{start} = $info2[2];
					 $result{$id}->{end} = $info2[3];
				         $result{$id}->{len} = abs($info2[3] - $info2[2])+1;
					 $result{$id}->{pattern_start} = $pattern[0];
		
					 } 
					 for my $keys (sort {	$result{$b}->{len} <=> $result{$a}->{len}
						        			 ||
								$result{$a}->{mis} <=> $result{$b}->{mis}
										 ||
								$result{$a}->{start}<=>$result{$b}->{start}
										||
								$result{$b}->{pattern_start} <=> $result{$a}->{pattern_start}

							    } keys %result){
					 my $len =$result{$keys}->{len};
					 my $temp_line = substr($line4,($result{$keys}->{start}-1),$len);
					 if((($result{$keys}->{pattern_len}==11)&&($len ==11)&&($result{$keys}->{mis} <=2)) || (($result{$keys}->{pattern_len}==10)&&($len ==10)&&($result{$keys}->{mis} <=2)) || (($result{$keys}->{pattern_len}==9)&&($len ==9)&&($result{$keys}->{mis} <=2)) || (($result{$keys}->{pattern_len}==8)&&($len ==8)&&($result{$keys}->{mis} <=1)) || (($result{$keys}->{pattern_len}==7)&&($len ==7)&&($result{$keys}->{mis} <=1))  ){
					print OUTPUT1 "$hash{$key}->{array}->[$i]->{name}\t$key\t$hash{$key}->{array}->[$i]->{start}:$hash{$key}->{array}->[$i]->{end}:$hash{$key}->{array}->[$i]->{str}:$hash{$key}->{array}->[$i]->{te_start}:$hash{$key}->{array}->[$i]->{te_end}\t$hash{$key}->{array}->[$j]->{start}:$hash{$key}->{array}->[$j]->{end}:$hash{$key}->{array}->[$j]->{str}:$hash{$key}->{array}->[$j]->{te_start}:$hash{$key}->{array}->[$j]->{te_end}\t$hash_seq{$result{$keys}->{id}}\t$temp_line\t$len\t$result{$keys}->{pattern_start}\t$result{$keys}->{start}\t$result{$keys}->{mis}\n"; 
					$mark{$i}=1;
					$mark{$j}=1;
					goto line1;
					 		}

						}		
				
#					}
					}else{
					goto line1;
					}
					
				}
					
			}
		
		}
		}
					line1: next;


	}
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
	  my $len = abs($end1 -$start2) +1;
	  if($len >=50){
	  $overlap=1;
	  }
          goto line2;
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
