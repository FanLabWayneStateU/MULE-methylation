#!/usr/bin/perl -w
use strict;
my $file = "MULE_gene/data/distance/rice_prob_chr_all";
open(INPUT,$file);
my $output = $file."_new";
open(OUTPUT,">$output");
while(my $line = <INPUT>){
chomp $line;
if($line =~ /^>/){
print OUTPUT "$line\n";
}else{
$line =~ s/\s+//g;
print OUTPUT "$line\n";
}
}
