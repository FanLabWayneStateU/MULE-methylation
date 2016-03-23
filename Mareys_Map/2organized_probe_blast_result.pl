#!/usr/bin/perl -w
use strict;
my $blast_out = "MULE_gene/data/distance/japo_prob_genome";
my $new_out = $blast_out."_organized";
`sort $blast_out -k 1,1 -k 11,11g $blast_out > $new_out`;
my $output = $blast_out."_single";
open(OUTPUT,">$output");
open(INPUT,$blast_out);
my $line_pre = <INPUT>;
chomp $line_pre;
my @info = split(/\s+/,$line_pre);
my $name_pre = $info[0];
my $e_val_pre = $info[10];
my $cdna_pre = substr($info[1],0,13);
my %already;
while(my $line = <INPUT>){
chomp $line;
@info = split(/\s+/,$line);
my $name = $info[0];
my $e_val = $info[10];
my $cdna = substr($info[1],0,13);
if(!(exists $already{$name})){
if(($name_pre ne $name) && ($e_val_pre <= 1e-10)){
print OUTPUT "$line_pre\n";
$already{$name_pre}=1;
}elsif(($name_pre eq $name) && ($e_val_pre <=1e-10)&&($e_val_pre < $e_val)){
print OUTPUT "$line_pre\n";
$already{$name_pre}=1;
}elsif(($name_pre eq $name)&&($e_val_pre <= 1e-10)&&($e_val_pre ==$e_val)&&($cdna_pre eq $cdna)){
$line_pre = $line;
$name_pre = $name;
$e_val_pre = $e_val;
$cdna_pre = $cdna;
}elsif(($name_pre eq $name)&&($e_val_pre <= 1e-10)&&($e_val_pre ==$e_val)&&($cdna_pre ne $cdna)){
$already{$name_pre}=1;
}elsif($e_val_pre > 1e-10){
$already{$name_pre}=1;
}
}
$name_pre = $name;
$e_val_pre = $e_val;
$cdna_pre = $cdna;
$line_pre = $line;
}


