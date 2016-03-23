#!/usr/bin/perl -w
use strict;
my $gff_file = "MULE_gene/data/japo_MULE_parental_seq_genome_blastn_gene_new1_exon30_parental";
my $expression_out = "$gff_file"."_gene_structure";

open(OUTPUT,">$expression_out");
open(INPUT,$gff_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $cds_name = $info[0];
my @coor1 = split(/:/,$info[0]);

for (my $i=1; $i <=$#info; $i++){
my @coor = split(/:/,$info[$i]);

my $cds_start = $coor[1];
my $cds_end = $coor[2];
my $cds_str = $coor[3];
my $chr_name1 = substr($coor[0],3);
if($chr_name1<10){
$chr_name1 = "0"."$chr_name1";
}
my $chr_name = "OsjapChr"."$chr_name1";


my ($exp_string,$exp_flag) = get_exp($chr_name,$cds_start,$cds_end,$cds_str);
my @exp_info = split(/\s+/,$exp_string);
my $new_exp_string = "";
if($coor1[1] eq $coor[0]){
for(my $i=0; $i<=$#exp_info; $i++){
my @id = split(/:/,$exp_info[$i]);
my $overlap = overlap($coor1[2],$coor1[3],$id[4],$id[5]);
if($overlap==0){
$new_exp_string .= "$exp_info[$i]\t";
}
}
}else{
$new_exp_string = $exp_string;
}

print OUTPUT "$cds_name\t$info[$i]\t$new_exp_string\n";
}

}


sub get_exp{

my ($chr_name,$cds_start,$cds_end,$cds_str)=(@_);
my $expression_file = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";
`awk '{if((\$1=="$chr_name")&&(\$3=="gene")){print }}' $expression_file | sort -k 4,4n -k 5,5n > gene_expression_glab4`;

open(INPUT1,"gene_expression_glab4");
my $exp_flag=0;
my $exp_str="";
my $mule_overlap=0;
my $exon_overlap=0;
my $rna_info;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_str = $info1[6];
my $rna_start = $info1[3];
my $rna_end = $info1[4];
my @rna_name=split(/;/,$info1[8]);
my $rna_name = $rna_name[0];
#if($rna_str eq $cds_str){
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		$exon_overlap = ($rna_end - $cds_start+1)/($rna_end - $rna_start+1);
		$mule_overlap = ($rna_end - $cds_start+1)/($cds_end - $cds_start+1);
 
		$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$cds_start:$rna_end";
		$exp_str .="$rna_info\t";
		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$exon_overlap =($cds_end - $cds_start + 1)/($rna_end - $rna_start + 1);
		$mule_overlap =1;
 		$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$cds_start:$cds_end";
		$exp_str .="$rna_info\t";

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		$mule_overlap = ($rna_end - $rna_start+1)/($cds_end - $cds_start +1);
		$exon_overlap = 1;
		$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$rna_end";
		$exp_str .="$rna_info\t";

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			$mule_overlap = ($cds_end - $rna_start+1)/($cds_end - $cds_start + 1);
			$exon_overlap = ($cds_end - $rna_start+1)/($rna_end - $rna_start + 1);

			$rna_info = "$rna_name:exon_$exon_overlap:mule_$mule_overlap:$rna_str:$rna_start:$cds_end";
			$exp_str .="$rna_info\t";

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}
	}
#}

}
line1: return($exp_str,$exp_flag);

}

sub overlap{

my ($cds_start,$cds_end,$rna_start,$rna_end)=(@_);
my $exp_flag=0;
my $exp_str="";
my $overlap;
my $rna_info;
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line2;

		}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		$overlap =1;
		goto line2;
		}elsif($rna_end <= $cds_start){
		$overlap=0;
		$exp_flag=0;
		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line2;


		}elsif(($rna_end >= $cds_end) &&($rna_start <= $cds_end)){
			$overlap = ($cds_end - $rna_start)/($cds_end - $cds_start + 1);
			$exp_flag=1;
			goto line2;

		}elsif(($rna_end > $cds_end) && ($rna_start > $cds_end)){
			$exp_flag=0;
			
		}
	}

line2: return($exp_flag);
}


