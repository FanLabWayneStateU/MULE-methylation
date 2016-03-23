#!/usr/bin/perl -w
use strict;
my $gene_gff = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff";


my $gene = "$gene_gff"."_gene";
`awk '{if(\$3 == "gene"){print}}' $gene_gff > $gene`;
my $output = $gene_gff."_non_TE_deep";

open(OUTPUT,">$output");
open(INPUT,"$gene");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $chr = substr($info[0],8);

if($chr =~ /(\d+)/){
$chr = $1;
if($chr =~ /^0(\d+)/){
$chr = $1;
}
}
$chr = "Chr$chr";
my $start = $info[3];
my $end = $info[4];

my $te = check_te($chr,$start,$end);
if($te ==0){
print OUTPUT "$line\n";
}
}

sub check_te{
my ($chr,$cds_start,$cds_end) = (@_);
my $te_gff = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new4.out.gff";
my $exp_flag=0;

`awk '{if(\$1=="$chr"){a=\$4-500; b=\$5+500; print a, b;}}' $te_gff  > gene_te_japo`;

open(INPUT1,"gene_te_japo");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
#my $rna_str = $info1[6];
my $rna_start = $info1[0];
my $rna_end = $info1[1];
#my $rna_name = $info1[9];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line1;
			}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		goto line1;

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line1;

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			goto line1;

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line2;
		  ;
		}
	}
}

line2:
my $te_gff1 = "MULE_gene/data/O.sativa.japonica_IRGSP_MSU_v7.fasta_new3.out.gff";
 

`awk '{if(\$1=="$chr"){a=\$4-1000; b=\$5+1000; print a,b ;}}' $te_gff1  > gene_te_japo`;

open(INPUT1,"gene_te_japo");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
#my $rna_str = $info1[6];
my $rna_start = $info1[0];
my $rna_end = $info1[1];
#my $rna_name = $info1[9];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line1;
			}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		goto line1;

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line1;

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			goto line1;

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line3;
		  ;
		}
	}
}

line3: 
my $te_gff2 = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5";
`awk '{if(\$2=="$chr"){split(\$3,a,":");split(\$4,b,":");c=a[1]-1000; d=b[2]+1000; print c,d;}}' $te_gff2 > gene_te_japo`;

open(INPUT1,"gene_te_japo");
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
#my $rna_str = $info1[6];
my $rna_start = $info1[0];
my $rna_end = $info1[1];
#my $rna_name = $info1[9];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exp_flag=1;
		goto line1;
			}elsif($rna_end >=$cds_end){
		$exp_flag=1;
		goto line1;

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$exp_flag=1;
		goto line1;

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$exp_flag=1;
			goto line1;

		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}

	}
}

line1: return($exp_flag);

}

