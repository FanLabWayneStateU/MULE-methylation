#!/usr/bin/perl -w
use strict;

my @target_spe = ("bart","brac","glum","meri","rufi","perr","japo","niva","glab","punc","indi");

for (my $k=0; $k<=0; $k++){

my $spe = $target_spe[$k];
my $temp_seq1 = "MULE_gene/data/$spe"."_para_temp1";
my $temp_seq2 = "MULE_gene/data/$spe"."_para_temp2";
my $temp_seq3 = "MULE_gene/data/$spe"."_para_temp3";
my $temp_seq4 = "MULE_gene/data/$spe"."_para_temp4";

my $mule_seq = "MULE_gene/data/$spe"."_MULE_internal_seq";
my %mule_seq = %{get_seq($mule_seq)};
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho";
my %ortho_record = %{get_ortho($file)};
my $coor ="MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
open(INPUTc,$coor);
my %coor;
while(my $linec = <INPUTc>){
chomp $linec;
my @infoc = split(/\s+/,$linec);
my @infoc1 = split(/:/,$infoc[0]);
my $id = "$infoc1[0]:$infoc1[1]:$infoc1[2]:$infoc1[3]";
$coor{$infoc[1]} = $id;
} 
my $list1 = "MULE_gene/data/$spe"."_MULE_para_1_list";
my $list2 = "MULE_gene/data/$spe"."_MULE_para_2_list";
my $list4 ="MULE_gene/data/$spe"."_MULE_para_4_list";
my $list6 ="MULE_gene/data/$spe"."_MULE_para_6_list";
my $list7 = "MULE_gene/data/$spe"."_MULE_para_7_list";
my $list8 = "MULE_gene/data/$spe"."_MULE_para_8_list";
my $list9 = "MULE_gene/data/$spe"."_MULE_para_9_list";
my $list1_dir = "MULE_gene/data/$spe"."_MULE_para_1_list_dir";
`mkdir $list1_dir`;
my $list2_dir = "MULE_gene/data/$spe"."_MULE_para_2_list_dir";
`mkdir $list2_dir`;
my $list4_dir = "MULE_gene/data/$spe"."_MULE_para_4_list_dir";
`mkdir $list4_dir`;
my $list6_dir = "MULE_gene/data/$spe"."_MULE_para_6_list_dir";
`mkdir $list6_dir`;
my $list7_dir = "MULE_gene/data/$spe"."_MULE_para_7_list_dir";
`mkdir $list7_dir`;
my $list8_dir = "MULE_gene/data/$spe"."_MULE_para_8_list_dir";
`mkdir $list8_dir`;
my $list9_dir = "MULE_gene/data/$spe"."_MULE_para_9_list_dir";
`mkdir $list9_dir`;
my $cur_dir;





my $psl_file = "MULE_gene/data/$spe"."_MULE_all_MULE.psl_best_sort";
open(INPUT,$psl_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $mule_id = $info[9];
my $ortho_line = $ortho_record{$coor{$mule_id}};
my @ortho = split(/\s+/,$ortho_line);
my $flag=0;
if(($ortho[1]=~/0/)&&($ortho[2]=~/0/)&&($ortho[3] =~ /0/)&&($ortho[4] =~ /0/)&&($ortho[5] =~ /0/)&&($ortho[6]=~/0/)&&($ortho[7]=~/0/)&&($ortho[8]=~/0/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list1");
$cur_dir = "$list1_dir/$mule_id";
open(OUTPUT1,">$cur_dir");
}elsif(($ortho[1]=~/0/)&&($ortho[2]=~/0/)&&($ortho[3] =~ /0/)&&($ortho[4] =~ /0/)&&($ortho[5] =~ /1/)&&($ortho[6]=~/0/)&&($ortho[7]=~/0/)&&($ortho[8]=~/0/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list2");
$cur_dir = "$list2_dir/$mule_id";
open(OUTPUT1,">$cur_dir");
}elsif(($ortho[1]=~/1/)&&($ortho[2]=~/1/)&&($ortho[3] =~ /1/)&&($ortho[4] =~ /1/)&&($ortho[5] =~ /1/)&&($ortho[6]=~/0/)&&($ortho[7]=~/0/)&&($ortho[8]=~/0/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list6");
$cur_dir = "$list6_dir/$mule_id";
}elsif(($ortho[1]=~/1/)&&($ortho[2]=~/1/)&&($ortho[3] =~ /1/)&&($ortho[4] =~ /1/)&&($ortho[5] =~ /1/)&&($ortho[6]=~/1/)&&($ortho[7]=~/0/)&&($ortho[8]=~/0/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list7");
$cur_dir = "$list7_dir/$mule_id";
}elsif(($ortho[1]=~/1/)&&($ortho[2]=~/1/)&&($ortho[3] =~ /1/)&&($ortho[4] =~ /1/)&&($ortho[5] =~ /1/)&&($ortho[6]=~/1/)&&($ortho[7]=~/1/)&&($ortho[8]=~/0/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list8");
$cur_dir = "$list8_dir/$mule_id";
}elsif(($ortho[1]=~/1/)&&($ortho[2]=~/1/)&&($ortho[3] =~ /1/)&&($ortho[4] =~ /1/)&&($ortho[5] =~ /1/)&&($ortho[6]=~/1/)&&($ortho[7]=~/1/)&&($ortho[8]=~/1/)&&($ortho[9]=~/0/)&&($ortho[10]=~/0/)){
open(OUTPUT,">>$list9");
$cur_dir = "$list9_dir/$mule_id";
}else{
$flag=1;
}
if($flag==0){
open(OUTPUT2,">$temp_seq1");
open(OUTPUT3,">$temp_seq2");
print OUTPUT2 ">query\n$mule_seq{$mule_id}\n";
print OUTPUT3 ">hit\n$mule_seq{$info[13]}\n";
open(OUTPUT1,">$cur_dir");
my ($mule_seq1,$target_seq1,$flag) = get_blat($temp_seq1,$temp_seq2,$line,$temp_seq3,$temp_seq4);
print OUTPUT1 ">$info[9]\n$mule_seq1\n";
print OUTPUT1 ">$info[13]\n$target_seq1\n";
print OUTPUT "$mule_id\n";
}
}
}

sub get_blat{
my ($temp_seq1,$temp_seq2,$line20,$temp_seq3,$temp_seq4) = (@_);

my @info20 = split(/\s+/,$line20);
my $str= $info20[8];
my $qlen = $info20[10];
my $qstart = $info20[11];
my $qend = $info20[12];
my ($tstart,$tend) = (0,0);
$tstart = $info20[15];
$tend = $info20[16];

my @chuck = split(/,/,$info20[18]);
my @mule_start = split(/,/,$info20[19]);
my @parent_start = split(/,/,$info20[20]);
my $temp_mule_seq="";
my $temp_parent_seq="";
if($str eq "+"){
for(my $j=0; $j <=$#mule_start; $j++){
my $temp_mule_start = $mule_start[$j]+1;
my $temp_mule_end = $mule_start[$j]+$chuck[$j]-1+1;
print "$temp_mule_start\t$temp_mule_end\n";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $temp_seq1 $temp_mule_start $temp_mule_end + > $temp_seq3`;
open(INPUT2,"$temp_seq3");
my $line2 = <INPUT2>;
chomp $line2;
$temp_mule_seq .= $line2;

}
for(my $j=0; $j <=$#parent_start; $j++){
my $temp_parent_start = $parent_start[$j]+1;
my $temp_parent_end = $parent_start[$j]+1+$chuck[$j]-1;
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $temp_seq2 $temp_parent_start $temp_parent_end + > $temp_seq4`;
open(INPUT3,"$temp_seq4");
my $line3 = <INPUT3>;
chomp $line3;
$temp_parent_seq .= $line3;
}
}elsif($str eq "-"){
for(my $j=$#mule_start; $j >=0; $j--){
my $temp_mule_end = $qlen-($mule_start[$j]);
my $temp_mule_start =  $qlen-($mule_start[$j]+$chuck[$j]-1);

`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $temp_seq1 $temp_mule_start $temp_mule_end + > $temp_seq3`;
open(INPUT2,"$temp_seq3");
my $line2 = <INPUT2>;
chomp $line2;
$temp_mule_seq .= $line2;
}


for(my $j=$#parent_start; $j >=0; $j--){
my $temp_parent_start = $parent_start[$j]+1;
my $temp_parent_end = $parent_start[$j]+1+$chuck[$j]-1;
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $temp_seq2 $temp_parent_start $temp_parent_end - > $temp_seq4`;
open(INPUT3,"$temp_seq4");
my $line3 = <INPUT3>;
chomp $line3;
$temp_parent_seq .= $line3;
}
}
line1: 
my $flag=0;
if(!($temp_mule_seq)){
;}else{
$flag=1;
}
return($temp_mule_seq,$temp_parent_seq,$flag);
}



sub get_ortho{
my $ortho_file = shift(@_);
my %ortho;
open(INPUTo,"$ortho_file");
while(my $lineo= <INPUTo>){
chomp $lineo;
my @info = split(/\s+/,$lineo);
$ortho{$info[0]} = $lineo;
}
return(\%ortho);
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
