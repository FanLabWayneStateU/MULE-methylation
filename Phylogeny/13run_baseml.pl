#!/usr/bin/perl -w
use strict;
my $spe = "$ARGV[0]";
#my @target_spe = ("1","2","4","6","7","8","9");
#my @target_spe = ("1","2","6","7","8","9");
#my @target_spe = ("1","7","8","9");
#my @target_spe = ("1","8","9");
my @target_spe = ("1","9");

my $output_dir = "MULE_gene/data/para_MULE_baseml_$spe" ;

`mkdir $output_dir`;
for (my $h=0; $h<=$#target_spe; $h++){
my $dir = "MULE_gene/data/$spe"."_MULE_para_"."$target_spe[$h]"."_list_dir";
my $spe_list = "MULE_gene/data/$spe"."_MULE_para_$target_spe[$h]"."_list";

my $output_list = "$spe_list"."_dist";
open(OUTPUT,">$output_list");
open(INPUT,$spe_list);
my $sub_out_dir = "$output_dir/$target_spe[$h]";
`mkdir $sub_out_dir`;
while(my $line = <INPUT>){
chomp $line;
my $seq_file = $dir."/$line";
my $tree_file = "MULE_gene/data/test_tree";
my $outfile ="$sub_out_dir/$line";  

my $ctl_text = qq{

      seqfile = $seq_file  *lysozymeSmall.txt
     treefile = $tree_file * lysozymeSmall.trees
      outfile = $outfile

        noisy = 2   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 2   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 4
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
	Mgene = 0
	clock = 0

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
	* ndata = 1
    fix_alpha = 0   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.5  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 5   * # of categories in the dG or AdG models of rates
	fix_rho = 1
	rho = 0. 
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time
       Small_Diff = 7e-6
       cleandata = 1
       * icode = 0
       * fix_blength = -1
};
my $ctl_file = "MULE_gene/data/test_$spe"."_4.ctl";

open(WRITEFILE,">$ctl_file");
print WRITEFILE "$ctl_text";
`paml4.7/bin/baseml $ctl_file`;
my $dist = get_dist($outfile);
if(!$dist){
print "$target_spe[$h]\t$line\n";
exit;
}
print OUTPUT "$line\t$dist\n";
}
}

sub get_dist{
my $file = shift(@_);
open(INPUT1,$file);
my $dist;
while(my $line1 = <INPUT1>){
chomp $line1;
print "enter\n";
if($line1 =~ /^Distances/){
<INPUT1>;
<INPUT1>;
<INPUT1>;
my $dist_line = <INPUT1>;
chomp $dist_line;
my @info1 = split(/\s+/,$dist_line);
my @info2 = split(/\(/,$info1[1]);
$dist = $info2[0];
}
}
return($dist);
}
