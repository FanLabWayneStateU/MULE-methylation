#!/usr/bin/perl -w
use strict;
my @spe = ("indi","bart","glab","rufi","niva","glum","meri","brac","punc","perr","japo","long");
for my $spe (@spe){
my $file = "MULE_gene/data/pack_MULE_candidate_$spe"."_jiang_repeatscout_new_tir_final3_name";
my $name = $file."_name";
`awk '{print \$1}' $file > $name`;
my $output = $file."_5";
`perl general_tool/check_num.pl $file 5 > $output`;
}
