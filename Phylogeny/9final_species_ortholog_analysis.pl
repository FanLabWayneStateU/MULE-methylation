#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Seq;

my @spe= ("rufi","indi","niva","bart","japo","glab","brac","glum","perr","punc","meri");
my @spe_gtf = ("oryza_rufipogon","oryza_sativa_indica","oryza_nivara","oryza_barthii","oryza_sativa_japonica","oryza_glaberrima","oryza_brachyantha","oryza_glumaepatula","leersia_perrieri","oryza_punctata","oryza_meridionalis");

my @spe_genome_fasta = ("Oryza_rufipogon_W1943_pseudo_v1.1","O.sativa.indica_9311","O.nivara_v1.0.pseudo","O.barthii_v1.3","O.sativa.japonica_IRGSP_MSU_v7","O.glaberrima.v1","O.brachyantha.v1.4","O.glumaepatula.v1.5","Lperr_v1.4","O.punctata_v1.2","O.meridionalis.v1.3");



for (my $i=10; $i<=10; $i++){
my $spe = $spe[$i];
my $total_output = "MULE_gene/data/pack_MULE_candidate_$spe"."_11_spe_ortho_species_genome_pep2";

open(OUTPUT,">$total_output");
my $genome_ortho = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_final5_new_no_transposase_e9_11_spe_ortho_species";

my $pep_all_file = "MULE_gene/data/pack_MULE_candidate_$spe[$i]"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_pep_seq";
my $focus_genome = "";
my %pep_hash = %{get_seq($pep_all_file)};
open(INPUT,"$genome_ortho");
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
print OUTPUT "$info[0]\t";
  for(my $j=1; $j<=10; $j++){
      my @spe_name = split(/:/,$info[$j]);
	my $genome_ortho = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final5_new_no_transposase_e9_11_spe_ortho_lineage_gene_list_flank_list"."_$spe_name[0]"."_ortho1";
      open(INPUT1,"$genome_ortho");
      my $chr_line = "MULE_gene/data/chr_line2";
      `grep "$info[0]" $genome_ortho | cat > $chr_line`;
      my ($genome_chr_flag,$chr,$start,$end) = check_chr("$chr_line");
      $genome_chr_flag=1;
      if($genome_chr_flag == 0){
      print OUTPUT "$spe_name[0]:0\t";
      }elsif($genome_chr_flag ==1){
	my $pep_ortho = "MULE_gene/data/pack_MULE_species_specific_$spe"."_$spe_name[0]"."_genome_pep2";
      my $pep_line = "MULE_gene/data/pep_line2";
      `grep "$info[0]" $pep_ortho | awk '{if(\$11<1e-10){print}}' | sort -k 11,11g | cat > $pep_line `;
      my $other_spe_gtf; 
      my $m=0;
      if($spe_name[0] =~ /rufi/){
      $m=0;
      }elsif($spe_name[0] =~ /indi/){
      $m=1;
      }elsif($spe_name[0] =~ /niva/){
      $m=2;
      }elsif($spe_name[0] =~ /bart/){
      $m=3;
      }elsif($spe_name[0] =~ /japo/){
      $m=4;
      }elsif($spe_name[0]=~ /glab/){
      $m=5;
      }elsif($spe_name[0] =~ /brac/){
      $m=6;
      }elsif($spe_name[0] =~ /glum/){
      $m=7;
      }elsif($spe_name[0] =~ /perr/){
      $m=8;
      }elsif($spe_name[0] =~ /punc/){
      $m=9;
      }elsif($spe_name[0] =~ /meri/){
      $m=10;
      } 
      
      $other_spe_gtf = "MULE_gene/data/maker_gff1/$spe_gtf[$m]".".maker.gff";
      my $other_spe_gtf1 = $other_spe_gtf."_gene";
      `grep "gene" $other_spe_gtf > $other_spe_gtf1`; 
      my $pep_coor_flag = check_pep("$pep_line",$chr,$start,$end,$other_spe_gtf);
       if($pep_coor_flag ==1){
       print OUTPUT "$spe_name[0]:1\t";
       }elsif($pep_coor_flag==0){
       my $focus_genome = "MULE_gene/data/$spe_genome_fasta[$m]".".fasta";
	my %focus_genome = %{get_seq($focus_genome)}; 
        my $chr_seq_new = $focus_genome{$chr};
        my $chr_test = "MULE_gene/data/genome_chr_test2";
	open(OUTPUT11,">$chr_test");
	print OUTPUT11 ">$chr\n";
	print OUTPUT11 "$chr_seq_new";
	my $start1 = $start - 200;
	my $end1  = $end + 200;
 	`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $chr_test $start1 $end1 + > MULE_gene/data/genome_seq_12_test2`;
	open(INPUT20,"MULE_gene/data/genome_seq_12_test2");
	my $line2= <INPUT20>;
	my $genome_seq_other = "MULE_gene/data/genome_seq_12_test2_new";
	open(OUTPUT20,">$genome_seq_other");
	print OUTPUT20 ">test\n";
	print OUTPUT20 "$line2";

       my $glimmer_test = "MULE_gene/data/other_MULE_Glimmer_test2";
       `/wsu/home/fe/fe49/fe4960/GlimmerHMM/sources/glimmerhmm $genome_seq_other -d GlimmerHMM/trained_dir/rice -o $glimmer_test`; 
       my ($pep_file,$pep_flag) = get_pep_seq($glimmer_test,$genome_seq_other);
       if($pep_flag==1){
       my $species_pep = "MULE_gene/data/species_pep2";
       open(OUTPUTsp,">$species_pep");
       print OUTPUTsp ">$info[0]\n";
       print OUTPUTsp "$pep_hash{$info[0]}\n";
       my $pep_len = length($pep_hash{$info[0]});
       my $pep_blast_out = "MULE_gene/data/pep_blast_out2";
       `./blast-2.2.20/bin/formatdb -i $pep_file -p T`;
   
       `./blast-2.2.20/bin/blastall -p blastp -d $pep_file -i $species_pep -o $pep_blast_out -m 8`; 
       my $hit_pep = check_de_novo($pep_blast_out,$pep_len);
       if($hit_pep ==1){
       print OUTPUT "$spe_name[0]:1\t";
       }elsif($hit_pep ==0){
       print OUTPUT "$spe_name[0]:0\t";
       }
       }elsif($pep_flag==0){
       print OUTPUT "$spe_name[0]:0\t";
       }
      }
     }
  }
print OUTPUT "\n";
}
}

sub check_de_novo{
my ($pep_blast_out,$pep_len) = (@_);
my $flag=0;
open(INPUTpb,"$pep_blast_out");
my %hash;
while(my $line1 = <INPUTpb>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
my $iden = $info1[2];
my $coverage = $info1[3]/$pep_len;
if($info1[10]< 1e-10){

$flag=1;
goto line10;
}
}

line10: return($flag);

}



sub get_pep_seq{

my ($glimmer_output,$genome_seq) = (@_);
open(INPUTgl,"$glimmer_output");
my $pep_seq_file = "MULE_gene/data/GlimmerHMM_pep_seq2";
open(OUTPUTgl,">$pep_seq_file");
for(my $i=1; $i<=9; $i++){
<INPUTgl>;
}
my $pep_flag=0;
my %hash;
my $line_name=1;
while(my $line = <INPUTgl>){
$pep_flag=1;
chomp $line;
if($line){
my @info = split(/\s+/,$line);
my $gene_name = "$line_name";

if(!(exists $hash{$gene_name}->{count})){
$hash{$gene_name}->{count}=0;
$hash{$gene_name}->{str}=$info[3];
}
$hash{$gene_name}->{array}->[$hash{$gene_name}->{count}]->{start} = $info[5];
$hash{$gene_name}->{array}->[$hash{$gene_name}->{count}]->{end} = $info[6];

$hash{$gene_name}->{count}++;

}else{

$line_name++;
}
}
my %cds_hash;
for my $key (sort keys %hash){
#print OUTPUTgl ">$key\n";
print "$key\t$hash{$key}->{count}\n";
if($hash{$key}->{str} eq "-"){
for(my $i=$hash{$key}->{count}-1; $i>=0; $i--){
my $start = $hash{$key}->{array}->[$i]->{start};
my $end = $hash{$key}->{array}->[$i]->{end};
print "$start\t$end\n";
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $genome_seq $start $end $hash{$key}->{str} > cds_seq_test_test2`;
open(INPUTcds,"cds_seq_test_test2");
my $linecds = <INPUTcds>;
chomp $linecds;
#print OUTPUTgl "$linecds";
$cds_hash{$key} .= "$linecds";
}
}elsif($hash{$key}->{str} eq "+"){
for(my $i=0; $i<=$hash{$key}->{count}-1; $i++){
my $start = $hash{$key}->{array}->[$i]->{start};
my $end = $hash{$key}->{array}->[$i]->{end};
`/wsu/apps/gnu-4.4.7/perl/perl-5.20.1/bin/perl epigene/scripts/38extract_seq.pl $genome_seq $start $end $hash{$key}->{str} > cds_seq_test_test2`;

open(INPUTcds,"cds_seq_test_test2");
my $linecds = <INPUTcds>;
chomp $linecds;
#print OUTPUTgl "$linecds";
$cds_hash{$key} .= "$linecds";
}

}
}

for my $key (sort keys %cds_hash){
my $pep_seq = get_pep($cds_hash{$key});
print OUTPUTgl ">$key\n";
print OUTPUTgl "$pep_seq\n";
}

return($pep_seq_file,$pep_flag);
}


sub get_pep{ 
my $new_line=shift(@_);
open(OUTPUT1,">temp_pep2");
print OUTPUT1 ">pep\n";
print OUTPUT1 "$new_line\n";
close(OUTPUT1);
my $in = Bio::SeqIO->new(-file => "temp_pep2",-format => 'Fasta');
my $chr_seq = $in->next_seq();
my $protein = $chr_seq->translate()->seq();
return($protein);
}




sub get_genome_seq{
my ($chr,$start,$end,$genome_seq)=(@_);
my $genome_seq_file = "MULE_gene/data/genome_seq_12_test2";
$start = $start - 100;
$end  = $end + 100;
`perl MULE_gene/scripts/83genomic_seq.pl $chr $start $end + $genome_seq > $genome_seq_file`;
return($genome_seq_file);
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


sub check_chr {
my $file = shift(@_);
open(INPUTc,"$file");
my $line=<INPUTc>;
chomp $line;
my @info = split(/\s+/,$line);
my @focus = split(/:/,$info[0]);
my $focus_chr = $focus[1];
if($focus_chr =~ /(\d+)/){
$focus_chr ="$1";
if($focus_chr =~ /0(\d+)/){
$focus_chr = "$1";
}
}
my @concern = split(/:/,$info[2]);
my $concern_chr = $concern[1];
if($concern_chr =~ /(\d+)/){
$concern_chr ="$1";
if($concern_chr =~ /0(\d+)/){
$concern_chr = "$1";
}
}
my $flag=0;
if($concern_chr == $focus_chr){
$flag=1;
}
return($flag,$concern[1],$concern[2],$concern[3]);
}
      
sub check_pep{
my ($pep_file,$blat_chr,$blat_start,$blat_end,$gene_gtf)=(@_);
my $pep_coor_flag=0;
if(open(INPUTp,$pep_file)){
for(my $i=0; $i<=1; $i++){
if(my $line = <INPUTp>){
chomp $line;
my @info = split(/\s+/,$line);
my @gene_name = split(/\./,$info[1]);
my $gene_name=$gene_name[0];
if(!($gene_name[0] =~ /perr/)){
$gene_name = "O".substr(lc($gene_name[0]),1);
}
`grep "$gene_name" $gene_gtf | cat > gene_coor2`;
open(INPUTgc,"gene_coor2") ;
my $linegc = <INPUTgc>;

chomp $linegc;
if(! $linegc){
print "$gene_name\n";
}
my @info1 = split(/\s+/,$linegc);

my $focus_chr = $blat_chr;
if($focus_chr =~ /(\d+)/){
$focus_chr ="$1";
if($focus_chr =~ /0(\d+)/){
$focus_chr = "$1";
}
}

my $focus_chr1 = $info1[0];
if($focus_chr1 =~ /(\d+)/){
$focus_chr1 ="$1";
if($focus_chr1 =~ /0(\d+)/){
$focus_chr1 = "$1";
}
}

if($focus_chr1 == $focus_chr){

$pep_coor_flag = overlap($blat_start,$blat_end,$info1[3],$info1[4]);
if($pep_coor_flag==1){
goto line1;
}
}
}
}
}
line1: return($pep_coor_flag);
}

sub overlap{
my ($start1,$end1,$start2,$end2)= (@_);
if($start1 > $end1){
my $temp =$start1;
$start1 = $end1;
$end1 = $temp;
}
if($start2 > $end2){
my $temp = $start2;
$start2 = $end2;
$end2 = $temp;
}
my $overlap=0;
if($end1 < $start2){
goto line2;
}elsif($end1 >= $start2){
	if($start1 < $start2){
	  if($end1 > $end2){
	  my $len = abs($end2 - $start2)+1;
	  if($len >=20){
	  $overlap=1;
	  }
	  goto line2;
	  }elsif($end1 <=$end2){
	  my $len = abs($end1 -$start2) +1;
	  if($len >=20){
	  $overlap=1;
	  }
	  goto line2;
	  }
        }elsif($start1 >= $start2){
            if($end1 < $end2){
	       my $len = abs($end1 - $start1)+1;
		if($len >= 20){
		$overlap=1;
		}
		goto line2;
	     }elsif($end1 >= $end2){
		if($start1 < $end2){
		  my $len = abs($end2 - $start1)+1;
		   if($len >=20){
		   $overlap=1;
                   }
                   goto line2;
                }elsif($start1 >= $end2){
		   goto line2;
		}
	      }
        }
}

line2: return($overlap);
}


       



