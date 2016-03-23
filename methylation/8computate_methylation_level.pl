#! user/bin/perl -w
use strict;
###
# perl this.pl coor_file p_value_file replicate_1 replicate_2 remove_file
# ####
my @file = ("$ARGV[0]");
#my $bi_file = "MULE_gene/data/japo_methy/O_sativa_japonica.qvalue.txt";
my $bi_file = "MULE_gene/data/all_methy/O_sativa_japonica/genome_matrix.O_sativa_japonica.q24.cov3.txt.meth";
my $context_file = "MULE_gene/data/japo_methy/SRR037416_23_se_pe_bedGraph";
my $flag = "$ARGV[1]";

for my $file(@file){
#my $file_methy = "$file"."$flag"."_methy5"; #for read=5
my $file_methy = "$file"."$flag"."_methy3_1"; #for read=3

open(OUTPUT,">$file_methy");
my $file_methy_stat = "$file_methy"."_stat";
open(OUTPUT1,">$file_methy_stat");
my %methy;
my %context;
#my $cover_file = "MULE_gene/data/japo_methy/genome_matrix.O_sativa_japonica.txt";
my $cover_file = "MULE_gene/data/all_methy/O_sativa_japonica/genome_matrix.O_sativa_japonica.txt";
open(INPUT5,$cover_file);
my %cover;
while(my $line5 = <INPUT5>){
chomp $line5;
my @info5 = split(/\s+/,$line5);
my $chr = "Chr$info5[0]";
my $post = "$info5[1]";
my @detail = split(/\//,$info5[4]);
$methy{$chr}->{$post}->{cov}=$info5[2];
$methy{$chr}->{$post}->{methy} = $detail[2];
$methy{$chr}->{$post}->{unmethy} = $detail[3];
}
open(INPUT4,$bi_file);
open(INPUT6,$context_file);
while(my $line4 = <INPUT4>){
chomp $line4;
my @info4 = split(/\s+/,$line4);
my $chr = "Chr$info4[0]";
my $post = "$info4[1]";
$methy{$chr}->{$post}->{qvalue}=$info4[3];

}
while(my $line6 = <INPUT6>){
chomp $line6;
my @info6 = split(/\s+/,$line6);
my $chr_id = substr($info6[0],3);
if($chr_id =~/^0(\d+)/){
$info6[0]  = "Chr$1";
}else{
$info6[0] = "Chr$chr_id";
}
$methy{$info6[0]}->{$info6[1]}->{context} = $info6[5];
}

print "done!\n";
open(INPUT,"$file");
my %pair;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my ($query_chr,$query_loc,$query_str,$query_bi,$query_context,$query_read);
$query_read=0;
my $pair = $info[0];
if($info[5] eq "C"){
$query_chr = "$info[1]";
$query_loc = $info[2];
$query_str = $info[3];
my $line1;
if(exists $methy{$query_chr}->{$query_loc}->{qvalue}){
if(($methy{$query_chr}->{$query_loc}->{methy}+$methy{$query_chr}->{$query_loc}->{unmethy})>=3){

 $query_bi = "me";
}else{
$query_bi = "unco";
}
}elsif((! (exists  $methy{$query_chr}->{$query_loc}->{qvalue}))&&(exists $methy{$query_chr}->{$query_loc}->{cov})){
if(($methy{$query_chr}->{$query_loc}->{methy}+$methy{$query_chr}->{$query_loc}->{unmethy})>=3){

 $query_bi = "unme";
}else{
$query_bi = "unco";
}
}elsif((! (exists  $methy{$query_chr}->{$query_loc}->{qvalue}))&&(!(exists $methy{$query_chr}->{$query_loc}->{cov}))){
$query_bi = "unco"
}

$query_context = $methy{$query_chr}->{$query_loc}->{context};


	if($query_bi eq "me"){
	print OUTPUT "$line\tme\t$query_context\n";
	$pair{$pair}->{me}->{$query_context}++;
	}elsif($query_bi eq "unme"){
	print OUTPUT "$line\tunme\t$query_context\n";
	$pair{$pair}->{unme}->{$query_context}++;
	}elsif($query_bi eq "unco"){
	print OUTPUT "$line\tunco\t$query_context\n";
	$pair{$pair}->{unco}->{$query_context}++;
        }

	}
}
my @context= ("CG","CHH","CHG");
for my $key (sort keys %pair){
for(my $i=0; $i <=2; $i++){
my $perc=0;
my ($methy, $unmethy,$unco)=(0,0,0);

if(exists $pair{$key}->{me}->{$context[$i]}){
$methy = $pair{$key}->{me}->{$context[$i]};
}
if(exists $pair{$key}->{unme}->{$context[$i]}){
$unmethy = $pair{$key}->{unme}->{$context[$i]};
}
if(exists $pair{$key}->{unco}->{$context[$i]}){
$unco = $pair{$key}->{unco}->{$context[$i]};
}
if(($methy + $unmethy) > 0){
$perc = $methy/($methy+$unmethy);
}
my $total=0;
my $covered =0;
if(($methy + $unmethy +$unco)>0){
$total = $methy + $unmethy + $unco;
$covered = ($methy + $unmethy)/$total;
}
print OUTPUT1 "$key\t$methy\t$unmethy\t$perc\t$covered\t$context[$i]\n";
}
}
}
