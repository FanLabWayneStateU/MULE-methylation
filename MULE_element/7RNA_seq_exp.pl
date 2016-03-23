#!/usr/bin/perl -w
use strict;
my @tissue = ("root","leaf","panicle");
my @spe = ("meri","brac","glab");
#my @spe = ("bart");
for my $spe (@spe){
#my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon_gff_exon_genic";
my $gff_file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor_involved_exon30_gff_CDS_genic";

my $mRNA = $gff_file."_mRNA";
`grep "mRNA" $gff_file > $mRNA`;
#my $gff_file = "MULE_gene/data/bart_MULE_parental_seq_genome_blastn_gene_new_parental";
my $expression_out = "$gff_file"."_gene_exp";
open(OUTPUT,">$expression_out");
open(INPUT,$mRNA);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @coor1 = split(/;/,$info[8]);
my $chr_name = $info[0];
my $cds_start = $info[3];
my $cds_end = $info[4];
my $genename = substr($coor1[0],3);

#my $genename = "$coor1[1]:$coor1[2]:$coor1[3]:$coor1[4]:$coor1[5]:$coor1[6]";
print OUTPUT "$genename\t$cds_start\t$cds_end\t";
for my $tissue(@tissue){
my ($exp,$exp_str) = get_exp($chr_name,$cds_start,$cds_end,$tissue,$spe);

print OUTPUT "$tissue:$exp:$exp_str\t";
}
print OUTPUT "\n";
}
}


sub get_exp{

my ($chr_name,$cds_start,$cds_end,$tissue,$spe)=(@_);
my $expression_file = "MULE_gene/data/cuff_exp/$spe"."_$tissue".".fpkm_tracking";
#print "$chr_name,$cds_start,$cds_end,$tissue,$spe\n"; 
`awk '{split(\$7,a,":");split(a[2],b,"-");if((a[1]=="$chr_name")&&((($cds_start<=b[1])&&($cds_end>=b[1]))||(($cds_start>=b[1])&&($cds_start<=b[2])))&&(\$13 == "OK")){print b[1],b[2],\$10}}' $expression_file | sort -k 1,1n -k 2,2n > gene_expression_cuff2`;

#exit;
open(INPUT1,"gene_expression_cuff2");
my $exp_flag=0;
my $exp_str=0;
my $exp_str1="";
my $mule_overlap=0;
my $exon_overlap=0;
my $rna_info;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $rna_start = $info1[0];
my $rna_end = $info1[1];
my $rna_exp = $info1[2];
	if($rna_start <= $cds_start){
		if(($rna_end <= $cds_end) &&($rna_end > $cds_start)){
		$exon_overlap = ($rna_end - $cds_start+1)/($rna_end - $rna_start+1);
		$mule_overlap = ($rna_end - $cds_start+1)/($cds_end - $cds_start+1);
		if($mule_overlap >=0.5){
		$exp_str1 .= "$rna_start:$rna_end:";
		$exp_str += $rna_exp*$mule_overlap;
		$exp_flag += $mule_overlap;
		} 
		}elsif($rna_end >=$cds_end){
		$exon_overlap =($cds_end - $cds_start + 1)/($rna_end - $rna_start + 1);
		$mule_overlap =1;
		$exp_str += $rna_exp*$mule_overlap;
		$exp_flag += $mule_overlap;
		$exp_str1 .= "$rna_start:$rna_end:";

		}
	}elsif($rna_start > $cds_start){
		if($rna_end <=$cds_end){
		$mule_overlap = ($rna_end - $rna_start+1)/($cds_end - $cds_start +1);
		$exon_overlap = 1;
		if($mule_overlap >=0.5){
		$exp_str += $rna_exp*$mule_overlap;
		$exp_flag +=$mule_overlap;
		$exp_str1 .= "$rna_start:$rna_end:";

		} 

		}elsif(($rna_end >= $cds_end) &&($rna_start < $cds_end)){
			$mule_overlap = ($cds_end - $rna_start+1)/($cds_end - $cds_start + 1);
			$exon_overlap = ($cds_end - $rna_start+1)/($rna_end - $rna_start + 1);
		if($mule_overlap >=0.5){
		$exp_str += $rna_exp*$mule_overlap;
		$exp_flag +=$mule_overlap;
		$exp_str1 .= "$rna_start:$rna_end:";
		} 
	
		}elsif(($rna_end > $cds_end) && ($rna_start >= $cds_end)){
			goto line1;
		  ;
		}
	}

}
line1: my $exp=0;
if($exp_flag>0){
$exp = $exp_str/$exp_flag;
}
return($exp,$exp_str1);

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



