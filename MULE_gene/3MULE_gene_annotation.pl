#!/usr/bin/perl -w
use strict;
my $spe_fa = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta";
my %seq = %{get_seq($spe_fa)};

my $file = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30";
#my $output = $file."_seq";
#open(OUTPUT,">$output");
my $output1 = $file."_gff_exon";
#my $output1 = $file."_gff_CDS";

open(OUTPUTgff,">$output1");
my $transpose = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9";

my %transpose =%{get_hash($transpose)};
my $auto_mule = $transpose."_auto_mule";
my %auto_mule = %{get_hash($auto_mule)};

sub get_hash{
my $infile = shift(@_);
open(INPUTi,$infile);
my %hash;
while(my $line = <INPUTi>){
my @info = split(/\s+/,$line);
my @start = split(/:/,$info[2]);
my @end = split(/:/,$info[3]);
my $name = "$info[0]:$info[1]:$start[0]:$end[1]:$start[2]";
$hash{$name}=1;
}
return(\%hash);
}
 
#my $output = $file."_CDS_seq";
#open(OUTPUT,">$output");
my $agi_annotation = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";

my $agi_annotation_mRNA = $agi_annotation."_CDS";
my $agi_annotation_exon = $agi_annotation."_exon";
`awk '{if(\$3=="CDS"){print}}' $agi_annotation > $agi_annotation_mRNA`;
my %CDS;
my %exon;
open(INPUT,"$agi_annotation_mRNA");
while(my $line = <INPUT>){
chomp $line;
my @info =split(/\s+/,$line);
my @name = split(/;/,$info[8]);
my $gene_name = substr($name[0],7);
if(exists $CDS{$gene_name}->{cds_num}){
$CDS{$gene_name}->{array}->[$CDS{$gene_name}->{cds_num}]->{start} = $info[3];
$CDS{$gene_name}->{array}->[$CDS{$gene_name}->{cds_num}]->{end} = $info[4];
$CDS{$gene_name}->{cds_num}++;
}else{
$CDS{$gene_name}->{cds_num}=0;
$CDS{$gene_name}->{strand}=$info[6];
$CDS{$gene_name}->{array}->[$CDS{$gene_name}->{cds_num}]->{start} = $info[3];
$CDS{$gene_name}->{array}->[$CDS{$gene_name}->{cds_num}]->{end} = $info[4];
$CDS{$gene_name}->{cds_num}++;
}
}

open(INPUT1,"$agi_annotation_exon");
while(my $line = <INPUT1>){
chomp $line;
my @info =split(/\s+/,$line);
my @name = split(/;/,$info[8]);
my $gene_name = substr($name[0],7);
if(exists $exon{$gene_name}->{cds_num}){
$exon{$gene_name}->{array}->[$exon{$gene_name}->{cds_num}]->{start} = $info[3];
$exon{$gene_name}->{array}->[$exon{$gene_name}->{cds_num}]->{end} = $info[4];
$exon{$gene_name}->{cds_num}++;
}else{
$exon{$gene_name}->{cds_num}=0;
$exon{$gene_name}->{strand}=$info[6];
$exon{$gene_name}->{array}->[$exon{$gene_name}->{cds_num}]->{start} = $info[3];
$exon{$gene_name}->{array}->[$exon{$gene_name}->{cds_num}]->{end} = $info[4];
$exon{$gene_name}->{cds_num}++;
}
}

my $Glimmer = "$transpose"."_Glimmer_sort";
my $Glimmer_mRNA = $Glimmer."_CDS";
`awk '{if(\$3 == "CDS"){print}}' $Glimmer > $Glimmer_mRNA`;
open(INPUT1,"$Glimmer_mRNA");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 =split(/\s+/,$line1);
my @name1 = split(/;/,$info1[8]);
my @coor = split(/:/,$name1[0]);
my $mule_name = "$coor[0]:$coor[1]:$coor[2]:$coor[3]:$coor[4]";
my $mule =  substr($mule_name,3);
 
my $gene_name1 = substr($name1[0],3,(length($name1[0])-7));
if(exists $CDS{$mule}->{$gene_name1}->{cds_num}){
$CDS{$mule}->{$gene_name1}->{array}->[$CDS{$mule}->{$gene_name1}->{cds_num}]->{start} = $info1[3];
$CDS{$mule}->{$gene_name1}->{array}->[$CDS{$mule}->{$gene_name1}->{cds_num}]->{end} = $info1[4];
$CDS{$mule}->{$gene_name1}->{cds_num}++;
}else{
$CDS{$mule}->{$gene_name1}->{cds_num}=0;
$CDS{$mule}->{$gene_name1}->{strand}=$info1[6];
$CDS{$mule}->{$gene_name1}->{array}->[$CDS{$mule}->{$gene_name1}->{cds_num}]->{start} = $info1[3];
$CDS{$mule}->{$gene_name1}->{array}->[$CDS{$mule}->{$gene_name1}->{cds_num}]->{end} = $info1[4];
$CDS{$mule}->{$gene_name1}->{cds_num}++;
}
}

open(INPUT2,"$file");
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
my @mule = split(/:/,$info2[0]);
my $chr = $mule[1];
#open(OUTPUTt,">temp_chr1_new");
#print OUTPUTt ">$chr\n";
#print OUTPUTt "$seq{$chr}\n";

if((exists $transpose{$info2[1]}) &&(! (exists $auto_mule{$info2[1]}))){
if(exists $info2[2]){
my %large;
    for(my $i=2; $i<=$#info2; $i++){
	my @info1 = split(/:/,$info2[$i]);
	my @gene = split(/\./,$info1[0]);
	my @mule_overlap = split(/_/,$info1[2]);

	if(exists $large{$gene[0]}){
		if($mule_overlap[1] > $large{$gene[0]}->{overlap}){
			$large{$gene[0]}->{line} = $info2[$i];
			$large{$gene[0]}->{overlap} = $mule_overlap[1];
		}
	}else{
		$large{$gene[0]}->{line} = $info2[$i];
		$large{$gene[0]}->{overlap} = $mule_overlap[1];
	}
    }
for my $key (sort keys %large){
my $new_line = $large{$key}->{line};
my @new_info = split(/:/,$new_line);
my $gene_name = $new_info[0];
my $chr = $mule[1];
#open(OUTPUTt,">temp_chr1_new");
#print OUTPUTt ">$chr\n";
#print OUTPUTt "$seq{$chr}\n";

my $str = $CDS{$gene_name}->{strand};
my @cds = @{$CDS{$gene_name}->{array}};
my @exon =@{$exon{$gene_name}->{array}}; 
my $seq_string="";
if($str eq "+"){
print OUTPUTgff "$mule[1]\tGlimmer_AGI\tmRNA\t$exon[0]->{start}\t$exon[$#exon]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name;Parent=$info2[1]:$gene_name\n";
}elsif($str eq "-"){
print OUTPUTgff "$mule[1]\tGlimmer_AGI\tmRNA\t$exon[$#exon]->{start}\t$exon[0]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name;Parent=$info2[1]:$gene_name\n";
}
for(my $j=0; $j<=$#cds; $j++){
#`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new $cds[$j]->{start} $cds[$j]->{end} $str > cds_seq1_new`;
#open(INPUTt,"cds_seq1_new");
#my $linet = <INPUTt>;
#chomp $linet;
#$seq_string .= "$linet";
}
#print OUTPUT ">$info2[1]:$gene_name\n";
#print OUTPUT "$seq_string\n";
for(my $j=0; $j<=$#exon; $j++){
#for(my $j=0; $j<=$#cds; $j++){

my $j_id = $j+1;
print OUTPUTgff "$mule[1]\tGlimmer_AGI\texon\t$exon[$j]->{start}\t$exon[$j]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name:exon.$j_id;Parent=$info2[1]:$gene_name\n";
#print OUTPUTgff "$mule[1]\tGlimmer_AGI\tCDS\t$cds[$j]->{start}\t$cds[$j]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name:cds.$j_id;Parent=$info2[1]:$gene_name\n";

}
}
}else{
my $mule_name = $info2[1];
for my $key (sort keys %{$CDS{$mule_name}}){
my $gene_name = $key;
$gene_name =~ s/cds/gene/g ;
my $str = $CDS{$mule_name}->{$key}->{strand};
my @cds = @{$CDS{$mule_name}->{$key}->{array}};
my $seq_string="";

if($str eq "+"){
#print OUTPUTgff "$mule[1]\tGlimmer_AGI\tmRNA\t$cds[0]->{start}\t$cds[$#cds]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name;Parent=$info2[1]:$gene_name\n";
print OUTPUTgff "$mule[1]\tGlimmer_AGI\tmRNA\t$cds[0]->{start}\t$cds[$#cds]->{end}\t.\t$str\t.\tID=$gene_name;Parent=$gene_name\n";

}elsif($str eq "-"){
#print OUTPUTgff  "$mule[1]\tGlimmer_AGI\tmRNA\t$cds[$#cds]->{start}\t$cds[0]->{end}\t.\t$str\t.\tID=$info2[1]:$gene_name;Parent=$info2[1]:$gene_name\n";
print OUTPUTgff  "$mule[1]\tGlimmer_AGI\tmRNA\t$cds[$#cds]->{start}\t$cds[0]->{end}\t.\t$str\t.\tID=$gene_name;Parent=$gene_name\n";

}


for(my $j=0; $j<=$#cds; $j++){
#`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl temp_chr1_new  $cds[$j]->{start} $cds[$j]->{end} $str > cds_seq1_new`;
#open(INPUTt,"cds_seq1_new");
#my $linet = <INPUTt>;
#chomp $linet;
#$seq_string .= "$linet";
my $j_id = $j+1;
print  OUTPUTgff "$mule[1]\tGlimmer_AGI\texon\t$cds[$j]->{start}\t$cds[$j]->{end}\t.\t$str\t.\tID=$gene_name:exon.$j_id;Parent=$gene_name\n";
#print  OUTPUTgff "$mule[1]\tGlimmer_AGI\tCDS\t$cds[$j]->{start}\t$cds[$j]->{end}\t.\t$str\t.\tID=$gene_name:cds.$j_id;Parent=$gene_name\n";

}

#print OUTPUT ">$gene_name\n";

#print OUTPUT "$seq_string\n";
}
}
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
