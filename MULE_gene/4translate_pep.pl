#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Seq;

my @spe = ("japo","bart","glum","indi","niva","rufi","meri","glab","punc","perr","brac");
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_CDS_seq";

my $pep = $file."_pep_screen";

open(OUTPUT,">$pep");
open(INPUT,$file);

my $id;
while(my $line = <INPUT>){
chomp $line;
if($line =~ /^>(\S+)/){
$id = "$1";
}elsif($line =~ /^[A-Z]/){
my $len = length($line);
if((($len % 3) ==0)&&($len >=150)){
my $protein = get_pep($line);
my $new_line = substr($protein,0,length($protein)-1);
if(($new_line =~ /\*/)||(!($new_line =~ /^M/))||(!($protein =~ /\*$/))){
;}else{
print OUTPUT ">$id\n";
print OUTPUT "$protein\n";
}

}
}
}
}
sub get_pep{ 
my $new_line=shift(@_);
open(OUTPUT1,">temp_pep");
print OUTPUT1 ">pep\n";
print OUTPUT1 "$new_line\n";
close(OUTPUT1);
my $in = Bio::SeqIO->new(-file => "temp_pep",-format => 'Fasta');
my $chr_seq = $in->next_seq();
my $protein = $chr_seq->translate()->seq();
return($protein);
}
