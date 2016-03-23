#!/usr/bin/perl -w
use strict;
my $coor_file = $ARGV[0];
my $coor_temp = $ARGV[1];
my $outputerr = $ARGV[2];
open(INPUTc,"$coor_file");
open(OUTPUTerr,">$outputerr");
my $dir = "MULE_gene/data/methy_slide_wd_f2";
`mkdir $dir`;
my $output;
if($coor_file =~ /data\/(\S+)$/){
$output = "$dir/$1";
print "$output\n";
}
my $file2 = "MULE_gene/data/japo_internal_MULE_c_coor_single_tir_fl_japo_methy_w_methy3_1_genic1";

open(INPUT2,$file2);
my %fl;
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
$fl{$info2[0]}->{$info2[1]}->{$info2[2]}->{context} = $info2[7];
$fl{$info2[0]}->{$info2[1]}->{$info2[2]}->{readnum} = $info2[8];
$fl{$info2[0]}->{$info2[1]}->{$info2[2]}->{methy} = $info2[6];
}
my $file = "MULE_gene/data/japo_internal_MULE_c_coor_single_tir_japo_methy_w_methy3_1_genic1";

open(INPUT,$file);
my %tir;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$tir{$info[0]}->{$info[1]}->{$info[2]}->{context} = $info[7];
$tir{$info[0]}->{$info[1]}->{$info[2]}->{readnum} = $info[8];
$tir{$info[0]}->{$info[1]}->{$info[2]}->{methy} = $info[6];
}
my $file1 = "MULE_gene/data/japo_internal_MULE_c_coor_single_new_japo_methy_w_methy3_1_genic1";

open(INPUT1,"$file1");
my %internal;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
$internal{$info1[0]}->{$info1[1]}->{$info1[2]}->{context} = $info1[7];
$internal{$info1[0]}->{$info1[1]}->{$info1[2]}->{readnum} = $info1[8];
$internal{$info1[0]}->{$info1[1]}->{$info1[2]}->{methy} = $info1[6];
}
my @TIR_L;
my @TIR_R;
my @internal_L;
my @internal_R;
my @FL_L;
my @FL_R;
open(OUTPUT,">$output");
my @context = ("CG","CHH","CHG");
my $line_num=0;
while(my $linec = <INPUTc>){
chomp $linec;
my @infoc = split(/\s+/,$linec);
my @mule = split(/:/,$infoc[0]);
print OUTPUTerr "$infoc[0]\n";
my ($mule,$chr,$start,$end)=($mule[0],$mule[1],$mule[2],$mule[3]);
my $mule_coor = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

$line_num++;
`grep "$mule" $mule_coor | grep "$chr" | grep "$start" | grep "$end" > $coor_temp`;
open(INPUT2,"$coor_temp");
my $line2 = <INPUT2>;
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @left = split(/:/,$info2[2]);
my @right = split(/:/,$info2[3]);
my $tir_name ="$info2[0]:$info2[1]:$left[1]:$right[0]:+"; 
my $mule_name = "$info2[0]:$info2[1]:$left[1]:$right[0]:+";
my %new_hash2;
if(exists $fl{$tir_name}->{$chr}){
 %new_hash2 = %{$fl{$tir_name}->{$chr}};
}
my %new_hash;
if(exists $tir{$tir_name}->{$chr}){
 %new_hash = %{$tir{$tir_name}->{$chr}};
}
my %new_hash1;
 if(exists $internal{$mule_name}->{$chr}){
 %new_hash1 = %{$internal{$mule_name}->{$chr}};
}
my $flag=0;
my $L_FL_index=0;
#for(my $start=$left[1]; ($start>=$left[0]) && ($flag==0); $start = $start - 40){
for(my $start=$left[0]-1; ($start>=($left[0]-500)) && ($flag==0); $start = $start - 25){

#	my $end = $start-80;
	my $end = $start-50;

	my %pair;
	   if($end <= ($left[0]-500)){
		$end = ($left[0]-500);
		$flag=1;
	    }
		for my $key (sort {$new_hash2{$b}<=>$new_hash2{$a}} keys %new_hash2){
			   if(($key <= $start) && ($key >= $end)){
				#	if($new_hash2{$key}->{readnum} >= 5){
#					if($new_hash2{$key}->{readnum} >= 3){
					if(($new_hash2{$key}->{methy} eq "me")||($new_hash2{$key}->{methy} eq "unme")){ 
					$pair{$new_hash2{$key}->{context}}->{$new_hash2{$key}->{methy}}++;
					}elsif($new_hash2{$key}->{methy} eq "unco"){
					$pair{$new_hash2{$key}->{context}}->{$new_hash2{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
			if(($total > 0 ) &&($covered >=0.5)){
		#	if(($total > 0 )){

			$FL_L[$L_FL_index]->{$context}->{me} += $methy;
			$FL_L[$L_FL_index]->{$context}->{unme} += $unmethy;
			$FL_L[$L_FL_index]->{$context}->{unco} += $unco;			
			$FL_L[$L_FL_index]->{$context}->{count}++;
		#	print OUTPUT "$tir_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tL_TIR\n";
			}			
		}
$L_FL_index++;
}







$flag=0;
my $L_TIR_index=0;
my $TIR_win_size = int((abs($left[1]-$left[0])+1)/5);
#for(my $start=$left[1]; ($start>=$left[0]) && ($flag==0); $start = $start - 40){
#for(my $start=$left[1]; ($start>=$left[0]) && ($flag==0); $start = $start - 25){
for(my $start=$left[1]; ($start>=$left[0]) && ($flag==0); $start = $start - $TIR_win_size){

#	my $end = $start-80;
#	my $end = $start-50;
	my $end = $start - $TIR_win_size*2;
	my %pair;
	   if($end <= $left[0]){
	   
		$end = $left[0];
		$flag=1;
	    }
		for my $key (sort {$new_hash{$b}<=>$new_hash{$a}} keys %new_hash){
			   if(($key <= $start) && ($key >= $end)){
#					if($new_hash{$key}->{readnum} >= 5){
#				if($new_hash{$key}->{readnum} >= 3){
					if(($new_hash{$key}->{methy} eq "me")||($new_hash{$key}->{methy} eq "unme")){ 

					$pair{$new_hash{$key}->{context}}->{$new_hash{$key}->{methy}}++;
					}elsif($new_hash{$key}->{methy} eq "unco"){
					$pair{$new_hash{$key}->{context}}->{$new_hash{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
			if(($total > 0 ) &&($covered >=0.5)){
#			if(($total > 0 )){

			$TIR_L[$L_TIR_index]->{$context}->{me} += $methy;
			$TIR_L[$L_TIR_index]->{$context}->{unme} += $unmethy;
			$TIR_L[$L_TIR_index]->{$context}->{unco} += $unco;			
			$TIR_L[$L_TIR_index]->{$context}->{count}++;
		#	print OUTPUT "$tir_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tL_TIR\n";
			}
			
		}
$L_TIR_index++;
}
$flag=0;
my $internal_index=0;
my $internal_win = int((abs($left[1]-$right[0])+1)/38);
#for(my $start=$left[1]; ($start<=$right[0]) && ($flag==0); $start = $start + 100){
#for(my $start=$left[1]; ($start<=$left[1]+499) && ($flag==0); $start = $start + 25){
for(my $start=$left[1]; ($start<=$right[0]) && ($flag==0); $start = $start + $internal_win){

#	my $end = $start+200;
#	my $end = $start+50;
	my $end = $start+ 2*$internal_win;
	my %pair;
	   if($end >= $right[0]){
		$end = $right[0];
		$flag=1;
	    }
		for my $key (sort {$new_hash1{$a}<=>$new_hash1{$b}} keys %new_hash1){
			   if(($key >= $start) && ($key <= $end)){
			#		if($new_hash1{$key}->{readnum} >= 5){
				#	if($new_hash1{$key}->{readnum} >= 3){
					if(($new_hash1{$key}->{methy} eq "me")||($new_hash1{$key}->{methy} eq "unme")){ 

					$pair{$new_hash1{$key}->{context}}->{$new_hash1{$key}->{methy}}++;
					}elsif($new_hash1{$key}->{methy} eq "unco"){
					$pair{$new_hash1{$key}->{context}}->{$new_hash1{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
			if(($total > 0 ) &&($covered >=0.5)){
#			if(($total > 0 )){

	#		print OUTPUT "$mule_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tinternal\n";
	#		}
			$internal_L[$internal_index]->{$context}->{me} += $methy;
			$internal_L[$internal_index]->{$context}->{unme} += $unmethy;
			$internal_L[$internal_index]->{$context}->{unco} += $unco;		
			$internal_L[$internal_index]->{$context}->{count}++;	
		#	}else{
#			print "$total\t$linec\n";
#			exit;
			}
		
		}
#print "$internal_index\n";
$internal_index++;
}

goto line1;
$flag=0;
my $internal_index=0;
#for(my $start=$left[1]; ($start<=$right[0]) && ($flag==0); $start = $start + 100){
for(my $start=$right[0]; ($start>=$right[0]-499) && ($flag==0); $start = $start - 25){

	my $end = $start-50;

	my %pair;
	   if($end <= ($right[0]-499)){
		$end = ($right[0]-499);
		$flag=1;
	    }



	for my $key (sort {$new_hash1{$a}<=>$new_hash1{$b}} keys %new_hash1){
			#   if(($key >= $start) && ($key <= $end)){
			   if(($key <= $start) && ($key >= $end)){

#					if($new_hash1{$key}->{readnum} >= 5){
				#	if($new_hash1{$key}->{readnum} >= 3){
					if(($new_hash1{$key}->{methy} eq "me")||($new_hash1{$key}->{methy} eq "unme")){ 

					$pair{$new_hash1{$key}->{context}}->{$new_hash1{$key}->{methy}}++;
					}elsif($new_hash1{$key}->{methy} eq "unco"){
					$pair{$new_hash1{$key}->{context}}->{$new_hash1{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
			if(($total > 0 ) &&($covered >=0.5)){
	#		print OUTPUT "$mule_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tinternal\n";
	#		}
			$internal_R[$internal_index]->{$context}->{me} += $methy;
			$internal_R[$internal_index]->{$context}->{unme} += $unmethy;
			$internal_R[$internal_index]->{$context}->{unco} += $unco;		
			$internal_R[$internal_index]->{$context}->{count}++;	
			}		
		}
#print "$internal_index\n";
$internal_index++;
}
line1:
$flag=0;
my $R_TIR_index=0;
my $win_size = int((abs($right[1]-$right[0])+1)/5);
#for(my $start=$right[0]; ($start<=$right[1]) && ($flag==0); $start = $start + 40){
#for(my $start=$right[0]; ($start<=$right[1]) && ($flag==0); $start = $start + 25){
for(my $start=$right[0]; ($start<=$right[1]) && ($flag==0); $start = $start + $win_size){

#	my $end = $start+80;
#	my $end = $start+50;
	my $end = $start + $win_size*2;
	my %pair;
	   if($end >= $right[1]){
		$end = $right[1];
		$flag=1;
	    }
		for my $key (sort {$new_hash{$a}<=>$new_hash{$b}} keys %new_hash){
			   if(($key >= $start) && ($key <= $end)){
#					if($new_hash{$key}->{readnum} >= 5){
#					if($new_hash{$key}->{readnum} >= 3){
					if(($new_hash{$key}->{methy} eq "me")||($new_hash{$key}->{methy} eq "unme")){ 

					$pair{$new_hash{$key}->{context}}->{$new_hash{$key}->{methy}}++;
					}elsif($new_hash{$key}->{methy} eq "unco"){
					$pair{$new_hash{$key}->{context}}->{$new_hash{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
		#	if(($total > 0 )){
		
			if(($total > 0 ) &&($covered >=0.5)){

		#	print OUTPUT "$tir_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tR_TIR\n";
		#	}
			$TIR_R[$R_TIR_index]->{$context}->{me} += $methy;
			$TIR_R[$R_TIR_index]->{$context}->{unme} += $unmethy;
			$TIR_R[$R_TIR_index]->{$context}->{unco} += $unco;			
			$TIR_R[$R_TIR_index]->{$context}->{count}++;
			}
		}
$R_TIR_index++;
}


$flag=0;
my $R_FL_index=0;
#for(my $start=$right[0]; ($start<=$right[1]) && ($flag==0); $start = $start + 40){
for(my $start=$right[1]+1; ($start<=($right[1]+500)) && ($flag==0); $start = $start + 25){

#	my $end = $start+80;
	my $end = $start+50;

	my %pair;
	   if($end >= ($right[1]+500)){
		$end = $right[1]+500;
		$flag=1;
	    }
		for my $key (sort {$new_hash2{$a}<=>$new_hash2{$b}} keys %new_hash2){
			   if(($key >= $start) && ($key <= $end)){
			#		if($new_hash2{$key}->{readnum} >= 5){
	#				if($new_hash2{$key}->{readnum} >= 3){
					if(($new_hash2{$key}->{methy} eq "me")||($new_hash2{$key}->{methy} eq "unme")){ 

					$pair{$new_hash2{$key}->{context}}->{$new_hash2{$key}->{methy}}++;
					}elsif($new_hash2{$key}->{methy} eq "unco"){
					$pair{$new_hash2{$key}->{context}}->{$new_hash2{$key}->{methy}}++;
					}

			    }
		}
		for my $context (@context){
			my ($methy,$unmethy,$unco,$perc,$total,$covered)=(0,0,0,0,0,0);
			if(exists $pair{$context}->{me}){
			$methy = $pair{$context}->{me};
			}
			if(exists $pair{$context}->{unme}){
			$unmethy = $pair{$context}->{unme};
			}
			if(exists $pair{$context}->{unco}){
			$unco = $pair{$context}->{unco};
			}
			if(($methy + $unmethy) > 0){
			$perc = $methy/($methy + $unmethy);
			}
			$total = $methy + $unmethy + $unco;
			if($total >0){
			$covered = ($methy + $unmethy)/$total;
			}
			
			if(($total > 0 ) &&($covered >=0.5)){
	#	if(($total > 0 )){

		#	print OUTPUT "$tir_name\t$chr\t$start\t$end\t$methy\t$unmethy\t$perc\t$covered\t$context\tR_TIR\n";
		#	}
			$FL_R[$R_FL_index]->{$context}->{me} += $methy;
			$FL_R[$R_FL_index]->{$context}->{unme} += $unmethy;
			$FL_R[$R_FL_index]->{$context}->{unco} += $unco;			
			$FL_R[$R_FL_index]->{$context}->{count}++;
			}
		}
$R_FL_index++;
}





}
my %k;
$k{CG}=0;
$k{CHG}=0;
$k{CHH}=0;
#%k = %{print_out(\@TIR_L,$output,"L_TIR",\%k,40,$line_num,1)};
#%k = %{print_out(\@internal,$output,"internal",\%k,100,$line_num,2)};
#%k = %{print_out(\@TIR_R,$output,"R_TIR",\%k,40,$line_num,2)};
##########
print "$line_num\n";
#%k = %{print_out(\@FL_L,$output,"L_FL",\%k,25,$line_num,1)};
 %k = %{print_out(\@FL_L,$output,"L_FL",\%k,1,$line_num,1)};

#%k = %{print_out(\@TIR_L,$output,"L_TIR",\%k,25,$line_num,1)};
%k = %{print_out(\@TIR_L,$output,"L_TIR",\%k,1,$line_num,1)};
%k = %{print_out(\@internal_L,$output,"internal",\%k,1,$line_num,2)};

#%k = %{print_out(\@internal_L,$output,"internal",\%k,25,$line_num,2)};
#%k = %{print_out(\@internal_R,$output,"internal",\%k,25,$line_num,1)};
%k = %{print_out(\@TIR_R,$output,"R_TIR",\%k,1,$line_num,2)};

#%k = %{print_out(\@TIR_R,$output,"R_TIR",\%k,25,$line_num,2)};
#%k = %{print_out(\@FL_R,$output,"R_FL",\%k,25,$line_num,2)};
%k = %{print_out(\@FL_R,$output,"R_FL",\%k,1,$line_num,2)};

sub print_out{

my ($array,$output,$id,$k,$sw,$line_num,$flag)=(@_);
#print "yes\n";
my %k = %{$k};
open(OUTPUT,">>$output");
my @TIR_L = @{$array};
my @context = ("CG","CHH","CHG");
if($flag==1){
for (my $i=$#TIR_L; $i>=0; $i--){
for my $context (@context){
		
my ($methy,$unmethy,$unco,$perc,$total,$covered,$nor_perc)=(0,0,0,0,0,0);

if(exists $TIR_L[$i]->{$context}->{me}){
$methy = $TIR_L[$i]->{$context}->{me};
}

if(exists $TIR_L[$i]->{$context}->{unme}){
$unmethy = $TIR_L[$i]->{$context}->{unme};
}

if(exists $TIR_L[$i]->{$context}->{unco}){
$unco = $TIR_L[$i]->{$context}->{unco};
}

if(($methy + $unmethy) > 0){
$perc = $methy/($methy + $unmethy);
#$nor_perc = $perc/$line_num;
}

$total = $methy + $unmethy + $unco;
if($total >0){
$covered = ($methy + $unmethy)/$total;
}
if($TIR_L[$i]->{$context}->{count}>0){
$nor_perc = $TIR_L[$i]->{$context}->{count}/$line_num;
}
#if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.05) ){

#if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.1)){
if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.2)){
#if(($total > 0 ) &&($covered >=0.5)){

print OUTPUT "$k{$context}\t$methy\t$unmethy\t$perc\t$covered\t$context\t$nor_perc\t$id\n";
$k{$context} += $sw;

#print "$k\t$methy\t$unmethy\t$perc\t$covered\t$context\t$id\n";
#exit;
}
}
}
}elsif($flag==2){
for (my $i=0; $i<=$#TIR_L; $i++){
for my $context (@context){
		
my ($methy,$unmethy,$unco,$perc,$total,$covered,$nor_perc)=(0,0,0,0,0,0);

if(exists $TIR_L[$i]->{$context}->{me}){
$methy = $TIR_L[$i]->{$context}->{me};
}

if(exists $TIR_L[$i]->{$context}->{unme}){
$unmethy = $TIR_L[$i]->{$context}->{unme};
}

if(exists $TIR_L[$i]->{$context}->{unco}){
$unco = $TIR_L[$i]->{$context}->{unco};
}

if(($methy + $unmethy) > 0){
$perc = $methy/($methy + $unmethy);
#$nor_perc = $perc/$line_num;
}

$total = $methy + $unmethy + $unco;
if($total >0){
$covered = ($methy + $unmethy)/$total;
}
if($TIR_L[$i]->{$context}->{count}>0){
$nor_perc = $TIR_L[$i]->{$context}->{count}/$line_num;
}
#if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.05)){
#if(($total > 0 ) &&($covered >=0.5)){

#if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.1)){
if(($total > 0 ) &&($covered >=0.5) && ($nor_perc >=0.2)){


print OUTPUT "$k{$context}\t$methy\t$unmethy\t$perc\t$covered\t$context\t$nor_perc\t$id\n";
#print "$k\t$methy\t$unmethy\t$perc\t$covered\t$context\t$id\n";
#exit;
$k{$context} += $sw;

}
}
}
}
return(\%k);
}

