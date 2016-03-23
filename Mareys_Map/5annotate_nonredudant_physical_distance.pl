#!/usr/bin/perl -w
use strict;

my $gff_file = "MULE_gene/data/maker_gff/oryza_japonica.maker.gff.txt_gene";
open(INPUT,$gff_file);
my %gff;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @gene_name = split(/[=|;]/,$info[8]);
$gff{$gene_name[1]}->{start} = $info[3];
$gff{$gene_name[1]}->{end} = $info[4];
$gff{$gene_name[1]}->{chrom} = $info[0];
}

my $mark_file= "MULE_gene/data/distance/chr_all_mark_multiple";
my %mark;
open(INPUT2,"$mark_file");
while(my $line2 = <INPUT2>){
chomp $line2;
$mark{$line2}=1;
}

    my $output1 = "MULE_gene/data/distance/chr_all"."_MareyMap_pure_nore";

    open(OUTPUT1,">$output1");
print OUTPUT1 "\"set\"\t\"map\"\t\"mkr\"\t\"phys\"\t\"gen\"\n";


for(my $i=1; $i<=12; $i++){
my %already;
    my $chr = "$i";
    if($i < 10){
    $chr = "0$i";
    }
   # my $file = "MULE_gene/data/distance/chr$chr"."_prob_id_gene_new";
  my $file = "MULE_gene/data/distance/chr$chr"."_prob_id_gene_new_pure";

    open(INPUT,$file);
#    my $output = "MULE_gene/data/distance/chr$chr"."_MareyMap";
    my $output = "MULE_gene/data/distance/chr$chr"."_MareyMap_pure";

    open(OUTPUT,">$output");
    print OUTPUT "\"set\"\t\"map\"\t\"mkr\"\t\"phys\"\t\"gen\"\n";
    while(my $line = <INPUT>){
    chomp $line;
    my @info = split(/\s+/,$line);
    my %hash;
    my ($mark,$gene_dist);
    if($line =~ /^F/){
    $mark = $info[3];
    $gene_dist = ($info[1]+$info[2])/2;
    for(my $j=4; $j<=$#info; $j++){
      my @gene = split(/\:/,$info[$j]);
      if(exists $gene[1]){
      my $gene = $gene[1];
      $hash{$gene}->{start} = $gff{$gene}->{start};
      $hash{$gene}->{end} = $gff{$gene}->{end};
      $hash{$gene}->{chrom} = $gff{$gene}->{chrom}
      }
    }
   }else{
   for(my $j=2; $j<=$#info; $j++){
     $mark = $info[1];
     $gene_dist = $info[0];
     my @gene = split(/\:/,$info[$j]);
     if(exists $gene[1]){
     my $gene = $gene[1];
     $hash{$gene}->{start} = $gff{$gene}->{start};
     $hash{$gene}->{end} = $gff{$gene}->{end};
     $hash{$gene}->{chrom} = $gff{$gene}->{chrom};
     }  
   }
   }
   if(!(exists $mark{$mark})){
    my ($start_pos,$end_pos,$chr);
    for my $key ( sort {$hash{$a}->{start} <=> $hash{$b}->{start}} keys %hash ){
      $start_pos = $hash{$key}->{start};
      $chr = $hash{$key}->{chrom};
      last;
    }

   for my $key ( sort {$hash{$b}->{end} <=> $hash{$a}->{end}} keys %hash ){
      $end_pos = $hash{$key}->{end};
      last;
    } 
   my $phy_dist = ($start_pos + $end_pos)/2;
    if(! exists $already{$chr}->{$phy_dist}->{$gene_dist}){
    print OUTPUT1 "\"Ojap\"\t\"$chr\"\t\"$mark\"\t$phy_dist\t$gene_dist\n";

    $already{$chr}->{$phy_dist}->{$gene_dist}=1;
    }
  }
}
}

           
