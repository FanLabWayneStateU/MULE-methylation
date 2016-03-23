#!/usr/bin/perl -w
use strict;
my @file = ("BCP","TCP","UNM");

for my $dir(@file){
my $file = "MULE_gene/data/sRNA_tophat_$dir/$dir".".bwa.sam";

my $output = "MULE_gene/data/sRNA_tophat_$dir/$dir"."_sort_sRNA_list_uniq_bwa0";
open(OUTPUT,">$output");

my $output1 = "MULE_gene/data/sRNA_tophat_$dir/$dir"."_sort_sRNA_list_two_more_bwa0";
open(OUTPUT1,">$output1");

my $count_file = "MULE_gene/data/sRNA/$dir"."_processed.txt_list";
my %count = %{get_count($count_file)};

open(INPUT,$file);
my %read;
while(my $line = <INPUT>){
chomp $line;
if($line =~ /^\d+/){
my @info = split(/\s+/,$line);

if($info[2] =~ /Chr/){
my $NM;
if($line =~ /NM\:i\:(\d+)/){
$NM="$1";
}
my $match_num =get_match($info[5]);

if($NM==0){
$read{$info[0]}->{count}++;
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{chrom}= $info[2];
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{start}= $info[3];
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{end} = $info[3] + $match_num-1;
$read{$info[0]}->{len} = $count{$info[0]}->{len};
$read{$info[0]}->{hit_num}= $count{$info[0]}->{count};
if(($info[11] eq "XT:A:U")){
;
}elsif($info[-1] =~ /Chr/){
my $uniq_num;
if($line =~ /X0\:i\:(\d+)/){
$uniq_num="$1";
}
my $new_line = substr($info[-1],5);
my @info1 = split(/;/,$new_line);
my $boundary;
if($#info1 < ($uniq_num -2)){
$boundary = $#info1;
}else{
$boundary = $uniq_num -2;
}
for(my $i=0; $i<=$boundary; $i++){
my @info2 = split(/,/,$info1[$i]);
if($info2[-1]==0){
$match_num = get_match($info2[-2]);
$read{$info[0]}->{count}++;
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{chrom}= $info2[0];
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{start}= substr($info2[1],1);
$read{$info[0]}->{array}->[$read{$info[0]}->{count}]->{end} = substr($info2[1],1) + $match_num -1;
}
}

}

}
}
}
}
for my $key (sort keys %read){
if($read{$key}->{hit_num}>0){
if($read{$key}->{count}==1){
print OUTPUT "$key\t$read{$key}->{array}->[1]->{chrom}\t$read{$key}->{array}->[1]->{start}\t$read{$key}->{array}->[1]->{end}\t$read{$key}->{hit_num}\t$read{$key}->{len}\n";
}elsif($read{$key}->{count}>=2){
for(my $i=1;$i<=$read{$key}->{count};$i++){
print OUTPUT1 "$key\t$read{$key}->{array}->[$i]->{chrom}\t$read{$key}->{array}->[$i]->{start}\t$read{$key}->{array}->[$i]->{end}\t$read{$key}->{hit_num}\t$read{$key}->{len}\n";
}

}
}
}
close(OUTPUT);
close(OUTPUT1);
}

sub get_match{
my $phrase = shift(@_);
my @info = split(/\D/,$phrase);
my @info1 = split(/\d+/,$phrase);
my $sum=0;
for(my $i=0; $i<=$#info; $i++){
if($info1[$i+1] ne "I"){
$sum += $1;
}
}
return($sum);
}

sub get_count{
my $file_name = shift(@_);
open(INPUTf,"$file_name");
my %hash;
while(my $linef = <INPUTf>){
chomp $linef;
my @info = split(/\s+/,$linef);
$hash{$info[0]}->{count}=$info[1];
$hash{$info[0]}->{len}=$info[2];

}
return(\%hash);
}
