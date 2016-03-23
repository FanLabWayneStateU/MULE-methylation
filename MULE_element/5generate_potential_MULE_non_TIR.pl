#!/usr/bin/perl -w
use strict;
my $seq_file = "MULE_gene/data/non_TIR_MULE.txt_sort";
my %seq_hash = %{get_seq($seq_file)};
my $file = "$seq_file".".out_non_mule";
my $output = "$seq_file"."_rm_other";
my $output1 = "miss_MULE";
open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
open(INPUTa,$file);
my %mule_agi;
while(my $linea = <INPUTa>){
chomp $linea;
my @info = split(/\s+/,$linea);
$mule_agi{$info[0]}=1;
}

for my $key (sort keys %seq_hash){
print $key;
my $subkey;
if(($key =~ /(\S+)L$/)||($key =~ /(\S+)R$/)){
$subkey = $1;
}
my $newkey = $subkey."L";
my $newkey1 = $subkey."R";
if((exists $mule_agi{$newkey}) || (exists $mule_agi{$newkey1})|| (exists $mule_agi{$key})){
;
}else{
print OUTPUT ">$key\n";
print OUTPUT "$seq_hash{$key}\n";


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
my $temp = $1;
$hash{$id} = $seq;
$id = $temp;
$seq ="";
}elsif($linefa =~ /^[A-Z|a-z]/){
$seq .= "$linefa";
}
}
$hash{$id} = $seq;
print "done!\n";
return(\%hash);
}

