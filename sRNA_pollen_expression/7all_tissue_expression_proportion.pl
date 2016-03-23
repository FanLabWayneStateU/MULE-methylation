#!/usr/bin/perl -w
use strict;
my %exp;
my @file = ("SRR074144","SRR074145","SRR074146","SRR074147","SRR074151","SRR074170","SRR074171");
for my $i (@file){
my $file = "MULE_gene/data/RNA_seq/$i/genes.fpkm_tracking";
$exp{$i}= get_exp($file);
}
my %new_exp;
for my $key (keys %{$exp{$file[0]}}){
	my $sum=0;
	for my $i (@file){
	$sum += $exp{$i}->{$key};
	}
	for my $i(@file){
	if($sum==0){
	$new_exp{$i}->{$key}=0;
	}else{
	$new_exp{$i}->{$key} =$exp{$i}->{$key}/$sum;

	}

	}
}
my %exp_note;
for my $i1 (@file){
my $file = "MULE_gene/data/RNA_seq/$i1"."_gff/genes.fpkm_tracking";
$exp_note{$i1}= get_exp($file);
}

my %new_exp_note;
for my $key (keys %{$exp_note{$file[0]}}){
	my $sum=0;
	for my $i (@file){
	$sum += $exp_note{$i}->{$key}
	}
	for my $i(@file){
	if($sum==0){
	$new_exp_note{$i}->{$key}=0;
	}else{
	$new_exp_note{$i}->{$key} =$exp_note{$i}->{$key}/$sum;

	}

	}
}


my $noTE = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_gene";

#my $noTE = "MULE_gene/data/maker_gff1/oryza_sativa_japonica.maker.gff_non_TE_deep_no_transposase_non_TE_internal";
my %noTE;
open(INPUT4,"$noTE");
while(my $line4 = <INPUT4>){
chomp $line4;
my @info = split(/\s+/,$line4);
my @info1 = split(/[=|;]/,$info[8]);
$noTE{$info1[1]}=1;
}



my $lineage = "MULE_gene/data/pack_MULE_candidate_japo_11_spe_ortho_species_genome_pep2";
`awk '{if((\$2~/0/)&&(\$3 ~/0/)&&(\$4 ~ /0/)&&(\$5~/0/)&&(\$6 ~/0/)&&(\$7 ~/0/)&&(\$8~/0/)&&(\$9~/0/)&&(\$10~/0/)&&(\$11~/0/)){print \$1}}' $lineage | cat > lineage_file`;
my $AA = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_AA";
my $Asian = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_Asian_all";

my $spe9 = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_BB_all";


my $coor = "MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_coor";
my %coor;
open(INPUT2,$coor);
while(my $line2 = <INPUT2>){
chomp $line2;
my @info2 = split(/\s+/,$line2);
$coor{$info2[1]} = substr($info2[0],0,(length($info2[0])-2));
}
$lineage="MULE_gene/data/pack_MULE_candidate_japo_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage"; 
my %lineage;
open(INPUT10,$lineage);
while(my $line10 = <INPUT10>){
chomp $line10;
my @info3 = split(/\s+/,$line10);
$lineage{$info3[0]}=1;
}

 
my %AA;
open(INPUT3,$AA);
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
$AA{$info3[0]}=1;
}


my %Asian;
open(INPUT3,$Asian);
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
$Asian{$info3[0]}=1;
}

my %spe9;
open(INPUT3,$spe9);
while(my $line3 = <INPUT3>){
chomp $line3;
my @info3 = split(/\s+/,$line3);
$spe9{$info3[0]}=1;
}

my %lineage_exp_ave;

for my $j (keys %exp){
($lineage_exp_ave{$j}->{ave},$lineage_exp_ave{$j}->{num})=get_ave1($new_exp{$j},\%coor,\%lineage);

}


my %all_exp_ave;
for my $h (keys %exp){
($all_exp_ave{$h}->{ave},$all_exp_ave{$h}->{num}) = get_ave2($new_exp{$h});
}

my %Asian_exp_ave;
for my $k (keys %exp){
($Asian_exp_ave{$k}->{ave},$Asian_exp_ave{$k}->{num}) = get_ave1($new_exp{$k},\%coor,\%Asian);

}

my %AA_exp_ave;
for my $k (keys %exp){
($AA_exp_ave{$k}->{ave},$AA_exp_ave{$k}->{num}) = get_ave1($new_exp{$k},\%coor,\%AA);


}


my %spe9_exp_ave;
for my $k (keys %exp){
($spe9_exp_ave{$k}->{ave},$spe9_exp_ave{$k}->{num}) = get_ave1($new_exp{$k},\%coor,\%spe9);
}

my %noTE_exp_ave;
for my $l (keys %exp_note){
($noTE_exp_ave{$l}->{ave},$noTE_exp_ave{$l}->{num}) = get_ave($new_exp_note{$l},\%noTE);
}

my $output = "MULE_gene/data/RNA_seq/japo_exp_all_new1";
open(OUTPUT,">$output");
print OUTPUT "\t\tlineage\tAsian\tAA\tBB\tall\tAll_gene\n";
for my  $i (@file){
print OUTPUT "$i\t$lineage_exp_ave{$i}->{ave}\t$Asian_exp_ave{$i}->{ave}\t$AA_exp_ave{$i}->{ave}\t$spe9_exp_ave{$i}->{ave}\t$all_exp_ave{$i}->{ave}\t$noTE_exp_ave{$i}->{ave}\n";
}

print OUTPUT "num\t$lineage_exp_ave{$file[0]}->{num}\t$Asian_exp_ave{$file[0]}->{num}\t$AA_exp_ave{$file[0]}->{num}\t$spe9_exp_ave{$file[0]}->{num}\t$all_exp_ave{$file[0]}->{num}\t$noTE_exp_ave{$file[0]}->{num}\n";
print OUTPUT "num\t$lineage_exp_ave{$file[1]}->{num}\t$Asian_exp_ave{$file[1]}->{num}\t$AA_exp_ave{$file[1]}->{num}\t$spe9_exp_ave{$file[1]}->{num}\t$all_exp_ave{$file[1]}->{num}\t$noTE_exp_ave{$file[1]}->{num}\n";
print OUTPUT "num\t$lineage_exp_ave{$file[2]}->{num}\t$Asian_exp_ave{$file[2]}->{num}\t$AA_exp_ave{$file[2]}->{num}\t$spe9_exp_ave{$file[2]}->{num}\t$all_exp_ave{$file[2]}->{num}\t$noTE_exp_ave{$file[2]}->{num}\n";

sub get_ave1{
my ($exp,$coor,$AA) = (@_);
my %exp = %{$exp};
my %coor = %{$coor};
my %AA = %{$AA};
my $sum=0;
my $n=0;

for my $key (keys %exp){
my @info = split(/:/,$key);
my $id = "$info[0]:$info[1]:$info[2]:$info[3]:$info[4]";

if(exists $AA{$coor{$id}}){
if($exp{$key}>0){
$sum++;
}
$n +=1;

}
}
my $ave = $sum/$n;
return($ave,$n);
}
sub get_hash{
my $file = shift(@_);
open(INPUT1,$file);
my %hash;
while(my $line1 = <INPUT1>){
chomp $line1;
$hash{$line1}=1;
}
return(\%hash);
}


sub get_ave2{
my ($hash) = shift(@_);
my %hash = %{$hash};
my $sum=0;
my $n=0;

for my $key(sort keys %hash){
if($hash{$key}>0){
$sum++;
}
$n +=1;

}
my $ave = $sum/$n;
return($ave,$n);

}




sub get_ave{
my ($hash,$hash1) = (@_);
my %lineage = %{$hash1};
my %hash = %{$hash};
my $sum=0;
my $n=0;

for my $key(sort keys %hash){
if(exists $lineage{$key}){
if($hash{$key}>0){
$sum++;
}
$n +=1;

}
}
my $ave = $sum/$n;
return($ave,$n);

}



sub get_exp{
my $file = shift(@_);
open(INPUT,$file);
<INPUT>;
my %exp;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$exp{$info[0]}=$info[9];
}
return(\%exp);
}
