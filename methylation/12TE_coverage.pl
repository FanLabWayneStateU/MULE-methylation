#!/usr/bin/perl -w
use strict;
my $spe = "japo";
#my $spe = "niva";
my $coor = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
open(INPUT3,"$coor");
my %coor;
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
#$coor{substr($info3[0],0,length($info3[0])-2)}=substr($info3[1],0,length($info3[1])-2);
$coor{substr($info3[1],0,length($info3[1])-2)}=substr($info3[0],0,length($info3[0])-2);

}
my $mask_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_internal_seq_only_screen.out.gff";
my $mask_sort = $mask_file."_sort";
`sort -k 1,1n -k 4,4n -k 5,5n $mask_file > $mask_sort`;
open(INPUT1,$mask_sort);
my %hash;
<INPUT1>;
<INPUT1>;
<INPUT1>;

while(my $line1 = <INPUT1>){
chomp $line1;
my @info = split(/\s+/,$line1);
my @name = split(/:/,$info[0]);
my $name = "$name[0]:$name[1]:$name[2]:$name[3]";
if(!(exists $hash{$name}->{count})){
$hash{$name}->{count}=0;
$hash{$name}->{array}->[$hash{$name}->{count}]->{start}= $info[3];
$hash{$name}->{array}->[$hash{$name}->{count}]->{end}= $info[4];
$hash{$name}->{count}++;
}else{
if(($info[3]>=$hash{$name}->{array}->[$hash{$name}->{count}-1]->{start}) &&($info[3]<=$hash{$name}->{array}->[$hash{$name}->{count}-1]->{end})){
if($info[4] > $hash{$name}->{array}->[$hash{$name}->{count}-1]->{end}){
$hash{$name}->{array}->[$hash{$name}->{count}-1]->{end} = $info[4];
}
}else{
 $hash{$name}->{array}->[$hash{$name}->{count}]->{start}= $info[3];
$hash{$name}->{array}->[$hash{$name}->{count}]->{end}= $info[4];
$hash{$name}->{count}++;
}
}
}
my $output = $mask_file."_cover";
open(OUTPUT,">$output");
for my $key (sort keys %hash){
my $sum=0;
my $line="";
for(my $i=0; $i<$hash{$key}->{count}; $i++){
$sum += abs($hash{$key}->{array}->[$i]->{end}- $hash{$key}->{array}->[$i]->{start})+1;
#$line .= "$hash{$key}->{array}->[$i]->{start}:$hash{$key}->{array}->[$i]->{end}"; 
}
my $internal = $coor{$key};
my @info = split(/:/,$internal);
my $length = $info[3]-$info[2]+1;
my $cover = $sum/$length;
print OUTPUT "$key\t$cover\t$line\n";
}
 
