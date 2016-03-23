#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;
#require "../../bin/zcj_fun.pm";
require "general_tool/zcj_fun.pm";
#################change###################choose one from the bellowing format
my $man = '';
my $help= '';
my $query_seq = '';
my $hit_seq = '';
my $output = '';
my $gff = '';
my $spe = '';
my $blat_dir = '';
my $blast2seq = '';
my $chrom = '';
my $codeml = '';
my $problem_loc = '';
my $detail ='';
my $rel = '';
my $blat_out = '';
my $kaks_file = '';
my $pair_list = '';
GetOptions('help|?' => \$help,
			man => \$man,
			'query_seq=s' => \$query_seq,
			'hit_seq=s' => \$hit_seq,
			'output=s' => \$output,
			'gff=s' =>\$gff,
			'spe=i' => \$spe,
			'blat=s' =>\$blat_dir,
			'blast2seq=s'=> \$blast2seq,
			'chrom=s' => \$chrom,
			'codeml=s' => \$codeml,
			'problem_loc=s' =>  \$problem_loc,
			'detail=s' => \$detail,
                        'rel=s' => \$rel,
			'blat_out=s' => \$blat_out,
			'kaks_file=s' => \$kaks_file,
			'pair_list=s' => \$pair_list
			)
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage(1) unless ($query_seq && $hit_seq && $gff &&  $spe && $output && $blat_dir && $blast2seq && $chrom && $codeml && $problem_loc && $detail && $rel && $blat_out && $kaks_file && $pair_list);

#my $blat_dir = $ARGV[0];
#my $blast2seq = $ARGV[1];
#my $unannotated_cds = $ARGV[2];
#my $annotated_cds = $ARGV[3];
#my $spe = $ARGV[4];
#my $gff = $ARGV[5];
#my $output = $ARGV[6];
#my $spe = $ARGV[7];
#my $chrom = $ARGV[8];
#######################################
#
#my %locs_seq1=&read_fasta_file("../Chr3.cds","LOC_Os03g01008.1","space");
#my %locs_seq2=&read_fasta_file("../../data-MSUv7/all.cds","LOC_Os01g01010.1","space");
my %locs_seq1 = %{get_seq($query_seq)};
my %locs_seq2 = %{get_seq($hit_seq)};
#my %jap_cds_id_info;
my %codon=&codon_table;
my @both_stop;
my @hit_stop;
my @query_stop;
my %predict_cds_seq;
my %query_cds_seq;
#################change#####################
my $jap_cds_id_info = get_gff($gff,$spe,$chrom);
my %jap_cds_id_info = %{$jap_cds_id_info};
###########################################
#open (GFF,"../../data-MSUv7/all.gff3");
#while (<GFF>){
#    next if $_!~/CDS/;
#    next if $_!~/Chr3/;
#    my @lines=split/\t/,$_;
#    my @tmp_id=split/\:/,$lines[8];
#    my @model_id=split/=/,$tmp_id[0];#ID=LOC_Os01g01010.2
#    my @cds_id=split/;/,$tmp_id[1];#cds_1;Parent=LOC_Os01g01010.2
#    $cds_id[0]=~s/cds_//;#cds_1
#    $jap_cds_id_info{$model_id[1]}{$cds_id[0]}=abs($lines[3]-$lines[4])+1;
#}

my @reverse;
my @frame_shift;
my $ll;
open (ERRORS,">$problem_loc") or die "wrong";
open (DET,">$detail") or die "wrong";
open (INDICA,">$rel") or die "wrong";
open(OUTPUTt,">$kaks_file") or die "wrong";

#my $blat_out = "blat_out1";
#my $blat_output = "test_blatout";
################change###############
#`$blat_dir $hit_seq $query_seq $blat_out`;
#my $new_blat_output = best_hit($blat_out);
#print "finish";
#exit;
######################################
#open (NEWKS,"blat_753_candidate_2glab_cds_find_para.query") or die "wrong";
#open (NEWKS,"arabidopsis/data/spe_mafft_trans_id") or die "wrong";
#open(NEWKS,"arabidopsis/data/ortho_absent_3_spe0_3spe_genome_u_u_trans_id") or die "wrong";
open(NEWKS,"$pair_list") or die "wrong";
#<NEWKS>;
#<NEWKS>;
#<NEWKS>;
#<NEWKS>;
#<NEWKS>;
while (<NEWKS>){
   # print "$_\n";
   chomp $_;
    my @lines=split/\s+/,$_;
  #  print "$lines[0]\t$lines[1]\n";
#    next if $lines[0]<100;
  #  next if $lines[1]/$lines[0]>0.10;
    my $loc1=$lines[0];
    my $loc2=$lines[1];
    next if $loc1 eq $loc2;
    my $model_id=$lines[0];
    my $model_id2=$lines[1];
    my $flag="";
    my $two_ids=$model_id."_".$model_id2;
    my @temp_key = keys %{$jap_cds_id_info{$model_id}};
    #print "@temp_key\n";
    
    #print "$two_ids\n";
    #exit;
    
    foreach my $cds (sort {$a<=>$b} keys %{$jap_cds_id_info{$model_id}}){
   #     print "inside:$cds\t$model_id\t$jap_cds_id_info{$model_id}\n";
	my $start_point=0 if $cds eq 1;
	my $tmp=$cds-1;
        if(($cds > 1) && (exists $jap_cds_id_info{$model_id}{$tmp})){
	$start_point+=$jap_cds_id_info{$model_id}{$tmp} ;
        }elsif(($cds > 1) && (!(exists $jap_cds_id_info{$model_id}{$tmp}))){
        $start_point = 0;
        }
       
        
	my $tmp_cds_seq=substr($locs_seq1{$model_id},$start_point,$jap_cds_id_info{$model_id}{$cds});
#	print "yes\n";
 #       print "$tmp_cds_seq\t$model_id\t$model_id2\t$cds\t$blast2seq\n";
  #      exit;
	$flag=&get_aligment($tmp_cds_seq,$locs_seq2{$model_id2},$model_id,$model_id2,$cds,$blast2seq);
        #print "$flag\n";
	#exit;
	last if $flag eq "both_rev";
    }

#    goto ERR if $flag eq "errors";
    goto ERR if $flag eq "both_rev";
    my $predict_gene_seq="";
    my $aligned_gene_seq="";
    foreach my $cds (sort {$a<=>$b} keys %{$predict_cds_seq{$two_ids}}){
	$predict_gene_seq.=$predict_cds_seq{$two_ids}{$cds};
	$aligned_gene_seq.=$query_cds_seq{$two_ids}{$cds};
    }
#    print "$two_ids\t$predict_gene_seq\t$aligned_gene_seq\n";
 #   exit;
    my ($all_tmp1,$all_tmp2,$all_aa1,$all_aa2,$yn_tmp1,$yn_tmp2)=&calculate_yn00($two_ids,$predict_gene_seq,$aligned_gene_seq);

    my $alignment_len=length($locs_seq1{$model_id});
    my $alignment_len1=length($locs_seq2{$model_id2});
    my $alignment_len2=length($yn_tmp1);
    print INDICA "$model_id\t$alignment_len1\t$model_id2\t$alignment_len2\t$alignment_len\n"; 
    print DET ">$two_ids\n$all_tmp1\n$all_tmp2\n$all_aa1\n$all_aa2\n";
    open (YN,">codmel1") or die "wrong";
    print  YN "\t2\t$alignment_len2\nquery\n$yn_tmp1\nhit\n$yn_tmp2\n";
    close YN;
###################
    open(CTL,">codmel1.ctl");
    my $ctl_text = qq{
    seqfile = codmel1
    treefile = codmel1.tree
    outfile = codmel1_rel
    noisy = 9
    verbose = 1
    runmode = -2
    seqtype = 1
    CodonFreq = 2 
    clock = 0
    model = 0
    NSsites = 0
    icode = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = 0
    omega = 0.1
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 10
    getSE = 0
    RateAncestor = 0
    method = 0
    Small_Diff = .5e-6
    };
    print CTL "$ctl_text";
    system "$codeml codmel1.ctl";
    unlink("2ML.dN");
    unlink("2ML.dS");
    unlink("2ML.t");
    unlink("2NG.dN");
    unlink("2NG.dS");
    unlink("2NG.t");
    system qq( tail codmel1_rel -n 1  > temp_codmel_output1);
    `tail codmel1_rel -n 4 | head -n 1 > temp_codmel_output2`;
    open(INPUTt,"temp_codmel_output1");
    my $line_codmel = <INPUTt>;
    chomp $line_codmel;
    open(INPUTt1,"temp_codmel_output2");
    my $line_codmel1 = <INPUTt1>;
    chomp $line_codmel1;
    my @infoa = split(/=/,$line_codmel1);
    my $ln_free = $infoa[1];
    #print "line:$line_codmel\n";
#    open(OUTPUTt,">>all_rels1.txt");
     open(OUTPUTt,">>$kaks_file");
    print OUTPUTt  "$two_ids\t$line_codmel\t$line_codmel1\n";
#exit;
###############################
 open(CTL,">codmel1.ctl");
    my $ctl_text = qq{
    seqfile = codmel1
    treefile = codmel1.tree
    outfile = codmel1_rel
    noisy = 9
    verbose = 1
    runmode = -2
    seqtype = 1
    CodonFreq = 2 
    clock = 0
    model = 0
    NSsites = 0
    icode = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = 1
    omega = 0.5
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 1
    getSE = 0
    RateAncestor = 0
    method = 0
    Small_Diff = .5e-6
    };
    print CTL "$ctl_text";
    system "$codeml codmel1.ctl";
    unlink("2ML.dN");
    unlink("2ML.dS");
    unlink("2ML.t");
    unlink("2NG.dN");
    unlink("2NG.dS");
    unlink("2NG.t");
    system qq( tail codmel1_rel -n 1  > temp_codmel_output1);
    `tail codmel1_rel -n 4 | head -n 1 > temp_codmel_output2`;
    open(INPUTt,"temp_codmel_output1");
    my $line_codmel = <INPUTt>;
    chomp $line_codmel;
    open(INPUTt1,"temp_codmel_output2");
    my $line_codmel1 = <INPUTt1>;
    chomp $line_codmel1;
    my @infoa = split(/=/,$line_codmel1);
    my $ln_05 = $infoa[1];
    my	$chi2 = 2*abs($ln_free - $ln_05);
    my 	$p = chi2_p($chi2,0,1);

    #print "line:$line_codmel\n";
#    open(OUTPUTt,">>all_rels1.txt");
     open(OUTPUTt,">>$kaks_file");
    print OUTPUTt  "$two_ids\t$line_codmel\t$line_codmel1\n";
    print OUTPUTt "$two_ids\tfixHalf:freePara\tchi2=$chi2\tdf=1\tp=$p\n";

#####################
 open(CTL,">codmel.ctl");
    my $ctl_text = qq{
    seqfile = codmel1
    treefile = codmel.tree
    outfile = codmel_rel
    noisy = 9
    verbose = 1
    runmode = -2
    seqtype = 1
    CodonFreq = 2 
    clock = 0
    model = 0
    NSsites = 0
    icode = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = 1
    omega = 1
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 1
    getSE = 0
    RateAncestor = 0
    method = 0
    Small_Diff = .5e-6
    };
    print CTL "$ctl_text";
    system "$codeml codmel.ctl";
    unlink("2ML.dN");
    unlink("2ML.dS");
    unlink("2ML.t");
    unlink("2NG.dN");
    unlink("2NG.dS");
    unlink("2NG.t");
    system qq( tail codmel_rel -n 1  > temp_codmel_output);
    `tail codmel_rel -n 4 | head -n 1 > temp_codmel_output1`;
    open(INPUTt,"temp_codmel_output");
    my $line_codmel = <INPUTt>;
    chomp $line_codmel;
    open(INPUTt1,"temp_codmel_output1");
    my $line_codmel1 = <INPUTt1>;
    chomp $line_codmel1;
    my @infoa = split(/=/,$line_codmel1);
    my $ln_1 = $infoa[1];
    my	$chi2 = 2*abs($ln_free - $ln_1);
    my 	$p = chi2_p($chi2,0,1);

    #print "line:$line_codmel\n";
#    open(OUTPUTt,">>all_rels1.txt");
     open(OUTPUTt,">>$kaks_file");
    print OUTPUTt  "$two_ids\t$line_codmel\t$line_codmel1\n";
    print OUTPUTt "$two_ids\tfixOne:freePara\tchi2=$chi2\tdf=1\tp=$p\n";
#exit;
ERR:
#exit if $two_ids eq "LOC_Os03g30990.1_LOC_Os01g01210.1";
}
print ERRORS "===reverse===\n".join("\n",@reverse);
print ERRORS "===shift===\n".join("\n",@frame_shift);
print ERRORS "===both-stop===\n".join("\n",@both_stop);
print ERRORS "===hit-stop===\n".join("\n",@hit_stop);
print ERRORS "===query-stop===\n".join("\n",@query_stop);
close ERRORS;
close INDICA;
close DET;
close NEWKS;

sub chi2_p {
open OUTFILE, ">tmp6.R";
my $values = "$_[0],$_[1]";
print OUTFILE
"data<-data.frame(x=c($values))\n",
"pvalues<-pchisq(data\$x,df=$_[2], lower.tail = FALSE)\n",
"write(pvalues,('tmp6.out'))";
close OUTFILE;
`/wsu/apps/gnu-4.4.7/R/R-3.1.1/bin/R < tmp6.R --no-save --slave`;
open INFILE, "<tmp6.out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm tmp6.R tmp6.out`;
close INFILE;
return($results[0]);
}

sub calculate_yn00(){
    my ($two_ids,$predict_seq,$aligned_seq)=@_;
    $predict_seq=~ tr/[a-z]/[A-Z]/;
    $aligned_seq=~ tr/[a-z]/[A-Z]/;
    my $yn_length=length($predict_seq);
    my $yn2_length=length($aligned_seq);
    if ($yn_length ne $yn2_length){
#	print "$predict_seq\n$aligned_seq\n";
#	print "lll length_not_ok\n";
#	print $two_ids;
#	exit;
 	print ERRORS "$two_ids\t===length-not-equal===\n";
	print ERRORS "$predict_seq\n$aligned_seq\n";
	print ERRORS "the alignment length is not equal to each other, please check\n";

   }
    my $new_predict_seq;
    my $new_aligned_seq;
    for (my $i=0;$i<$yn2_length;$i++){
	my $a=substr($aligned_seq,$i,1);
	my $b=substr($predict_seq,$i,1);
	if ($a eq "+"){
	    my $new_a=substr($aligned_seq,$i,3);
	    $a="";
	    $b="";
	    if ($new_a eq "+++"){
		$i=$i+2;
	    }else{
		push (@frame_shift,$two_ids);
	    }
	}
	$new_predict_seq.=$b;
	$new_aligned_seq.=$a;
    }
    if (substr($new_aligned_seq,-3)!~/-/){
	if ($codon{substr($new_aligned_seq,-3)} eq "*"){
	    $new_predict_seq=substr($new_predict_seq,0,-3);
	    $new_aligned_seq=substr($new_aligned_seq,0,-3);
	}
    }
    $yn_length=length($new_predict_seq);
    print ERRORS "$two_ids can't devided by 3, wrong length\n" if ($yn_length/3 ne int ($yn_length/3));#7839 b) LOC_Os03g27140.12511LOC_Os01g14570.119705526 c) LOC_Os03g01080.12511LOC_Os01g14570.118805403
    my $yn_tmp1;
    my $yn_tmp2;
    my $all_tmp1;
    my $all_tmp2;
    my $all_aa1;
    my $all_aa2;
    for (my $i=0;$i< $yn_length;$i=$i+3){
	my $tmp1=substr($new_predict_seq,$i,3);
	my $tmp2=substr($new_aligned_seq,$i,3);#
	last if length($tmp1) ne "3" and length($tmp2) ne "3";#
	if (($tmp1 =~/\-|\~|n|N/ and $tmp1!~/\+/) or $tmp2 =~/\-|\~|N/ ){
	    next;
	}elsif($tmp1 =~/\+/){
	    my $shift_a=0;
	    for (my $j=0;$j<3;$j++){
		my $shift=substr($tmp1,$j,1);
		$shift_a++ if $shift eq "+";
	    }
	    if ($shift_a eq 3){
		next;
	    }else{
		push (@frame_shift,$two_ids);
		next;
	    }
	}
	if ($codon{$tmp1} eq "*" or $codon{$tmp2} eq "*"){#
	    if ($codon{$tmp1} eq "*" and $codon{$tmp2} eq "*"){
		push (@both_stop,$two_ids);
	    }elsif($codon{$tmp1} eq "*" and $codon{$tmp2} ne "*"){
		push (@hit_stop,$two_ids);
	    }elsif($codon{$tmp1} ne "*" and $codon{$tmp2} eq "*"){
		push (@query_stop,$two_ids);
	    }
	}
	if ($codon{$tmp1} ne "*" and $codon{$tmp2} ne "*"){#
	    $yn_tmp1.=$tmp1;
	    $yn_tmp2.=$tmp2;
	}
	$all_tmp1.=$tmp1;
	$all_tmp2.=$tmp2;
	$all_aa1.="-".$codon{$tmp1}."-";
	$all_aa2.="-".$codon{$tmp2}."-";#
    }
    return ($all_tmp1,$all_tmp2,$all_aa1,$all_aa2,$yn_tmp1,$yn_tmp2);
}
my %multiple;
my %chimeric;

sub get_aligment()
{
    my ($dna_seq1,$dna_seq2,$mid1,$mid2,$cds,$blast2seq)=@_;
#    my $temp_seq2 = $dna_seq2;
 #   	    open(OUTPUTtemp,">output_temp");
#	    print OUTPUTtemp ">test\n";
#	    print OUTPUTtemp "$dna_seq2\n";
#	    close OUTPUTtemp;

    open (TMP,">tmpdna1.fa") or die "wrong";
    print TMP ">$mid1"."_query\n$dna_seq1\n";
    close TMP;
    open (TMP,">tmpdna2.fa") or die "wrong";
    print TMP ">$mid2"."_hit\n$dna_seq2\n";
    close TMP;
#exit;
    system "rm tmp.fas.txt" if (-e "tmp.fas.txt");
#    system "bl2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
##################change######################
    system "$blast2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -F F -o tmp.fas.txt";
###########################################
    my $split_hit=0;
    my $two_ids=$mid1."_".$mid2;
    my @multiple_identity;
    undef @multiple_identity;    
    my $max_identity=0;
    my $max_hit=0;
    my $max_match_length=0;
    my $max_len_hit=0;
    my $max_match_bit=0;
    my $max_len_bit=0;
    my $max_bit=0;
    my $newrev_flag=0;
NEWREV:
    open (TMP,"tmp.fas.txt") or die "wrong";#
    while(<TMP>){
	if ($_=~/_hit/){
	    my @lines=split/\t/,$_;
	    $split_hit++ if $lines[10] < 4e-5;
	    next if $lines[10] > 4e-5;
	    $multiple_identity[$split_hit-1]=$lines[2] if $split_hit>0;
	    if ($split_hit eq 1){
		$max_identity =$lines[2];
		$max_match_length=$lines[3];
		$max_hit=$split_hit;
		$max_len_hit=$split_hit;
		$max_match_bit=$lines[11];
		$max_len_bit=$lines[11];
	    }else{
		if  ($max_identity <$lines[2]){
		    $max_identity =$lines[2];
		    $max_hit=$split_hit;
		    $max_match_bit=$lines[11];
		}
		if ($max_match_length<$lines[3]){
		    $max_match_length=$lines[3];
		    $max_len_hit=$split_hit;
		    $max_len_bit=$lines[11];
		}
	    }
	}
    }
    my $follow="F";
    if ($split_hit>1){#more than one hit
	print ERRORS "more than one hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
	foreach my $tmp_line2 (@multiple_identity){
	    if
 (abs($max_identity - $tmp_line2)>5){#前，我们只关心identity最高的
		$chimeric{$two_ids}{$cds}++;
	    }
	}
	if (abs($max_match_bit-$max_len_bit)<50){
 if (exists $chimeric{$two_ids}{$cds}){
		$max_bit=$max_match_bit;
	    }else{
		$max_bit= $max_len_bit;
	    }
	}else{
	    $max_bit=$max_match_bit;
	    if ($max_match_bit<$max_len_bit){
		$max_bit=$max_len_bit;
		
	    }
	    $follow="T";
	}

    }

    close TMP; 
    open (TMP,"tmp.fas.txt") or die "wrong";
    my $new_tmp_split=0;
    while(<TMP>){
	if ($_=~/_hit/){
	    my @lines=split/\t/,$_;
	    $new_tmp_split++ if $lines[10] < 4e-5;
	    next if $lines[10] > 4e-5;
	    if ($follow eq "T"){
		next if $lines[11] ne "$max_bit";
	    }else{
		next if (exists $chimeric{$two_ids}{$cds} and $new_tmp_split ne $max_hit);#的区间
		next if (not exists $chimeric{$two_ids}{$cds} and $new_tmp_split ne $max_len_hit);#
	    }
	    if ($lines[6]<$lines[7] and $lines[9]<$lines[8]){#环
		my $rcvalue = reverse $locs_seq2{$mid2};
                $rcvalue =~ tr/ACGTacgt/TGCAtgca/;  
		open (TMP,">tmpdna2.fa") or die "wrong";
		print TMP ">$mid2"."_hit\n$rcvalue\n";
		close TMP;
	#	system "bl2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
#################change#################################
			system "$blast2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -F F -o tmp.fas.txt";
######################################################
		$newrev_flag++;
		if ($newrev_flag eq "1"){
		    goto NEWREV;
		}else{
		    print ERRORS "both direction have hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
		    return "both_rev";
		}
	    }
	}
    }
    close TMP;

    if ($split_hit eq "0"){
	print ERRORS "no good hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
	for (my $i=1;$i<=length($dna_seq1);$i++){  
	    $predict_cds_seq{$two_ids}{$cds}.="-";
	    $query_cds_seq{$two_ids}{$cds}.="-";
	}
	return "not_hit";
    }
    my $rev_flag="F";

REV:
    my $new_split=0;
    open (TMP,"tmp.fas.txt") or die "wrong";
    while(<TMP>){
	if ($_=~/_hit/){
	    my @lines=split/\t/,$_;
	    $new_split++ if $lines[10] < 4e-5;
	    next if $lines[10] > 4e-5;
	    next if $new_split>1;
	    $max_bit=$lines[11];
	    if ($lines[6]<$lines[7] and $lines[9]<$lines[8]){#死循环
		my $rcvalue = reverse $locs_seq2{$mid2};
                $rcvalue =~ tr/ACGTacgt/TGCAtgca/;  
		open (TMP,">tmpdna2.fa") or die "wrong";
		print TMP ">$mid2"."_hit\n$rcvalue\n";
		close TMP;
#		system "bl2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
################change##############################		
		system "$blast2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -F F -o tmp.fas.txt";
####################################################
		$locs_seq2{$mid2}=$rcvalue;

		push (@reverse,$two_ids);
		if ($rev_flag eq "F"){
		    $rev_flag = "T";
		    goto REV;#向而重新计算
		}else{
		    print ERRORS "both direction have hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
		    return "both_rev";
		}
	    }
    	    if ($lines[5] >0){#存在gap
		my $read_line=int(5*((($lines[7]-$lines[6])/60)+1)+11)+1;
		my $read_line2=int(5*((($lines[9]-$lines[8])/60)+1)+11)+1;
		$read_line=$read_line2 if $read_line<$read_line2;

		&have_gap($dna_seq1,$dna_seq2,$read_line,$cds,$two_ids,$lines[6],$max_bit,$blast2seq);#
		return "gap";
	    }
#    	    my $value=substr($locs_seq2{$mid2},($lines[8]-1),($lines[9]-$lines[8]+1));#
	    open(INPUTtemp,"tmpdna2.fa");
            <INPUTtemp>;
	    my $temp_seq2 = <INPUTtemp>;
	    my $value=substr($temp_seq2,($lines[8]-1),($lines[9]-$lines[8]+1));#
#    	    my $value=substr($locs_seq2{$mid2},21608000,5000);#

	    my $temp_start = $lines[8]-1;
	    my $temp_end = $lines[9] - $lines[8] +1;
#	    print "value1:$temp_start:$temp_end:$value\n";
#	    exit;
	    for (my $i=1;$i<$lines[6];$i++){
		$value="-".$value;#的开头部分，用”-“代替
	    }
	    while (length($value)<length($dna_seq1)){#
		$value.="-";
	    }
	    $predict_cds_seq{$two_ids}{$cds}=$value;
#	    print "value:$value\n";
	    $query_cds_seq{$two_ids}{$cds}=$dna_seq1;
#	    print "dna_seq1:$dna_seq1\n";
	}
    }
    return "ok";
}


sub have_gap()
{
    my ($dna_seq1,$dna_seq2,$read_line,$cds,$two_ids,$dash,$max_bit,$blast2seq)=@_;
    chomp $max_bit;
    while($max_bit=~s/ //){}
#    system "bl2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 0 -o tmp.fas.gap.txt";
###################change##################
     system "$blast2seq -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 0 -F F -o tmp.fas.gap.txt";
############################################   
    open (TMP,"tmp.fas.gap.txt") or die "wrong";
    my $count_line=0;
    my $query_seq;
    my $hit_seq;
    my $al_flag="F";
    while(<TMP>){
	if ($max_bit ne "0"){
	    if ($_=~/Score/){
		my @arr=split/=/,$_;
		while ($arr[1]=~s/  / /){}
		my @arrb=split/ /,$arr[1];
		if ($arrb[1] eq $max_bit){
		    $read_line-=8;
		    $al_flag="T";
		}else{
		    next;
		}
	    }
	    next if $al_flag eq "F";
	}
	$count_line++;
	last if $count_line>$read_line;
	while($_=~s/  / /){}
	my @lines=split/ /;
	if ($_=~/Sbjct/){
	    $hit_seq.=$lines[2];
	}
	if ($_=~/Query:/){
	    $query_seq.=$lines[2];
	}
    }
    my $tmp_seq=$query_seq;
    my $tmp_dash=0;
    while ($tmp_seq=~s/-//){
	$tmp_dash++;
    }
    $hit_seq=~s/-/+/g;
    $query_seq=~s/-/+/g;
    for (my $i=1;$i<$dash;$i++){#$dash=$lines[6]
	$hit_seq="-".$hit_seq;
	$query_seq="-".$query_seq;
    }
    while(length($query_seq)<(length($dna_seq1)+$tmp_dash)){

	$hit_seq.="-";
	$query_seq.="-";
    }
    $predict_cds_seq{$two_ids}{$cds}=$hit_seq;
    $query_cds_seq{$two_ids}{$cds}=$query_seq;
}
###################change#######################
sub get_gff{
my ($gff_file,$spe,$chrom)=(@_);
my %jap_cds_id_info;
open(GFF,$gff);
if (($spe ==1) || ($spe==2) || ($spe ==3)){ #gtf file 1:human;2:mouse; 3:fruitfly
while (<GFF>){
    next if $_!~/CDS/;
    my @lines=split/\s+/,$_;
	    if($lines[0] eq $chrom){	
	    my $tmp_id=substr($lines[11],1,(length($lines[11])-3));
    	    my $cds_id=substr($lines[13],1,(length($lines[13])-3));#cds_1;Parent=LOC_Os01g01010.2
	    #if($tmp_id eq "FBtr0078185"){
            #print "$tmp_id\t$cds_id\n";
            #}
	    #exit;
            $jap_cds_id_info{$tmp_id}{$cds_id}=abs($lines[3]-$lines[4])+1;
            }
}


}elsif($spe==4){ #gff file rice;
my %tmp_cds_hash;
while (<GFF>){
    chomp $_;
    next if $_!~/CDS/;
    #next if $_!~/$chrom/;
    my @lines=split/\t/,$_;
    my @tmp_id=split/\:/,$lines[8];
    #print $tmp_id[0];
    #exit;
    my @model_id=split/=/,$tmp_id[0];#ID=LOC_Os01g01010.2
    $tmp_cds_hash{$model_id[1]}++;
#    my @cds_id=split/;/,$tmp_id[1];#cds_1;Parent=LOC_Os01g01010.2
 my   $cds_id= $tmp_cds_hash{$model_id[1]};#cds_1
    #print "$model_id[1]\t$cds_id\n";
    #exit;
    $jap_cds_id_info{$model_id[1]}{$cds_id}=abs($lines[3]-$lines[4])+1;
}
}elsif($spe==5){ #gff file of arabidopsis; Chr1	phytozome8_0	CDS	3760	3913	.	+	0	ID=PAC:19656964.CDS.1;Parent=PAC:19656964;pacid=19656964
my %tmp_cds_hash;
while (<GFF>){
    chomp $_;
    next if $_!~/CDS/;
#    next if $_!~/$chrom/;
    my @lines = split(/\s+/,$_);
#    my @tmp_id = split(/[:|\.|;/,$lines[8]);
    #my @tmp_id = split(/\./,$lines[8]);
     my @tmp_id = split(/;/,$lines[8]);
     my @model_id = split(/[=|:]/,$tmp_id[0]);
#     my $model_id = substr($tmp_id[0],3,(length($tmp_id[0])-3));
     my $model_id = "$model_id[1]:$model_id[2]:$model_id[3]:$model_id[4]:$model_id[5]:$model_id[6]";
#     print $model_id;
#     exit;
#     my $cds_id	= $tmp_id[1];
     if(!(exists $tmp_cds_hash{$model_id})){
     $tmp_cds_hash{$model_id}=1;
     }
     my $cds_id = $tmp_cds_hash{$model_id};
     $jap_cds_id_info{$model_id}{$cds_id}=abs($lines[3]-$lines[4])+1;
     $tmp_cds_hash{$model_id}++;
}

}
return(\%jap_cds_id_info);
}
########
sub best_hit{
my $blat_out = shift(@_);
my $new_blat_out = "new_blatout";
open(OUTPUTb,">$new_blat_out");
open(INPUTb,$blat_out);
my %hash;
for(my $i=0; $i<=4; $i++){
my $temp_line = <INPUTb>;
chomp $temp_line;
print  OUTPUTb "$temp_line\n";
}
my %count;
while(my $line = <INPUTb>){
	chomp $line;
	my @info = split(/\s+/,$line);
	print "line:$line\n";
	my ($qcover,$tcover,$iden) = pslCal($line);
	my $iden_score = 100.0 - $iden * 0.1;
	my $blat_score = blat_score($line);
	if(!(exists $count{$info[9]})){
	$count{$info[9]}=0;
	}
	$hash{$info[9]}->{$count{$info[9]}}->{iden_score} = $iden_score;
	$hash{$info[9]}->{$count{$info[9]}}->{score} = $blat_score;
	$hash{$info[9]}->{$count{$info[9]}}->{qcover} = $qcover;
	$hash{$info[9]}->{$count{$info[9]}}->{tcover} = $tcover;
	$hash{$info[9]}->{$count{$info[9]}}->{line} = $line;
	$count{$info[9]}++;
}
my %count_key;
for my $key (sort keys %hash){
	my %sub_hash = %{$hash{$key}};
#	sub by_score {
#	$sub_hash{$b}->{score} <=> $sub_hash{$a}->{score}
#	||
#	$sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score}
#	}
	for my $sub_key (sort {	$sub_hash{$b}->{score} <=> $sub_hash{$a}->{score}
	||
	$sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score} } keys %sub_hash){
		if(!(exists $count_key{$key})){
		print  OUTPUTb "$sub_hash{$sub_key}->{line}\n";
		$count_key{$key} =1;
		}
	}
}
close OUTPUTb;
#open(INPUTtest,"new_blatout");
#while(my $linet = <INPUTtest>){
#print "$linet";
#}
return("new_blatout");
}

sub blat_score{
my $psl = shift(@_);
my @info = split(/\s+/,$psl);
#my $score = ($info[0] + $info[2]/2) - $info[1] - $info[5] - $info[7];
my $score = ($info[0] + $info[2]/2) - $info[1] - $info[4] - $info[6];

return($score);
}

sub pslCal{
my $psl = shift(@_);
my $sizeMul = 1;
my ($qAliSize,$tAliSize,$aliSize);
my $milliBad=0;
my $sizeDif;
my $insertFactor;
my $total;
my @info = split(/\s+/,$psl);

$qAliSize = $sizeMul * ($info[12] - $info[11]);
$tAliSize = $info[16] - $info[15];

my $qcover = $qAliSize/$info[10];
my $tcover = $tAliSize/$info[14];
if($qAliSize <= $tAliSize){
$aliSize = $qAliSize;
}else{
$aliSize = $tAliSize;
}

if($aliSize <= 0){
	goto line1;
}
$sizeDif = $qAliSize - $tAliSize;
	if($sizeDif < 0){
#		$sizeDif = -$sizeDif;
	$sizeDif=0;
       }
#$insertFactor = $info[5];
$insertFactor = $info[4];

#	$insertFactor += $info[7];
$total = ($sizeMul * ($info[0] + $info[1]+$info[2]));
if($total !=0){
$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + ceil(3*log(1+$sizeDif))))/$total;
line1: return ($qcover,$tcover,$milliBad);
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
######

__END__

=head1 NAME

new_method_caculate_kaks.pl - align DNA seq to compute ka/ks

=head1 SYNOPSIS

new_method_caculate_kaks.pl [options]

 Options:
   -help            brief help message
   -man             full documentation
   -query_seq       query cds sequence file
   -hit_seq         hit sequence file
   -output          output file
   -gff             for human, mouse, fruitfly using ensembl gtf format, for rice and arabidopsis using gff3 format
   -spe             spe: 1 human 2 mouse 3 fruitfly 4 rice 5 arabidopsis 
   -blat	    blat directory
   -blast2seq	    blast2seq directory
   -codeml	    codeml of PAML directory
   -chrom	    chrom of the query sequence 
   -problem_loc     the file records the problem locus
   -detail          the file record the DNA alignment and Protein alignment
   -rel 
   -blat_out        the output file of blat

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-query_seq>


>seq_name1

ATGCGGG...

=item B<-hit_seq>

>seq_name2

ATGGCT...

=item B<-output>

Output data to a file.

=item B<-gff>

indicate the annotation file of query sequence. for human, mouse, fruitfly use gtf file, for rice and arabidopsis use gff file

=item B<-spe>

indicate species index, 1:human 2:mouse 3:fruitfly 4:rice 5:arabidopsis

=item B<-codeml>

indicate the directory of codeml

=item B<-blat>

indicate the directory of blat

=item B<-blast2seq_dir>

indicate the directory of blast2seq

=item B<-chrom>

indicate the chromosome of query sequences.

=item B<-problem_loc>

the file recording the locus with problem

=item B<-detail>

the file recording the CDS and peptide alignment

=item B<-rel>
the file recording the rel information

=item B<-blat_out>

the output file of blat

=item B<-blat_out>

the file recording the kaks value

=back

=head1 DESCRIPTION

B<new_method_caculate_kaks> will ...

=head1 AUTHOR 

Chengjun Zhang @ University of Chicago
Jun	Wang	@ Wayne state university


=head1 BUGS

none.

=head1 COPYRIGHT 

This program is free software. You may copy or redistribute it under the same terms as Perl itself.

=cut

################################################

