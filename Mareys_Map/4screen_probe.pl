#!/usr/bin/perl -w
use strict;
for(my $i=1; $i<=12; $i++){
my $chr_id="$i";
   if($i<10){
	$chr_id = "0$i";
    }
my $file = "MULE_gene/data/distance/chr$chr_id"."_prob_id_gene";
my $output = $file."_new_pure";
open(INPUT,$file);
open(OUTPUT,">$output");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $flag=0;
my $gene_chr="00";

if($info[0] =~ /^F/){
   for(my $j=4; $j<=$#info; $j++){
       if($info[$j] =~ /Osjap(\d+)g/){
         $gene_chr = "$1";
         if($gene_chr ne $chr_id){
		$flag++;
         }
        }
    }
}else{
   for(my $j=2; $j<=$#info; $j++){
       if($info[$j] =~ /Osjap(\d+)g/){
         $gene_chr = "$1";
         if($gene_chr ne $chr_id){
		$flag++;
         }
        }
    }
}
   if($gene_chr eq "00"){
   $flag++;
   }
if($flag==0){
print OUTPUT "$line\n";
}
}
}

