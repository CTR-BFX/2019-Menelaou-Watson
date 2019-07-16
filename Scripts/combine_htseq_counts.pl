#!/usr/bin/perl


# Program name:        combine_htseq_counts_v01.pl
# Program version:     v01
# Author:              Malwina Prater

# combine_htseq_counts_v01.pl <ARGV[0]> <ARGV[1]> >  ARGV[2]
# MP_combine_abundances_kallisto.pl <abundance.txt> <abundance.txt>  <OUT>.txt

# perl MP_combine_abundances_kallisto_v01.pl SLX-9615.A010.C66Y8ANXX.s_5.r_1_kallisto_output/abundance.tsv SLX-9612.A023.C6H2YANXX.s_8.r_1_kallisto_output/abundance.tsv SLX-9615.A016.C66Y8ANXX.s_5.r_1_kallisto_output/abundance.tsv /Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0002/write_perl_script_abundances_combined/output_test_file.tsv
# perl MP_combine_abundances_kallisto_v01.pl 

##### write function for all bash:
##### for i in *count1.txt; do echo $i; perl combine_htseq_counts_v01.pl $i ${i/count1.txt/count2.txt} > ${i/count1.txt/combined_counts.txt} ; done



use strict;
use warnings;

use File::Basename;




my %hash_ENST;
my $line_number=0;

open(FILE1, $ARGV[0]) || die "Can't open file1 for reading: $!\n";
while(<FILE1>)
 {
    chomp();
    my $line = $_;
    
    my ($target_id, $est_counts) = split(/\t/, $line) ;

    if( $line_number >= 0)   # \s == wide space; \s+ == any space; [0-9]+ == any integer, any times; ([0-9]+) == put in special brackets so pts in into special variable $1. if (()) then $2.
      {
        $hash_ENST{$target_id} = $est_counts;
      } 

    $line_number++;
 }
close(FILE1);


open(FILE2, $ARGV[1]) || die "Can't open file2 for reading: $!\n";
$line_number=0;
while(<FILE2>)
 {
    chomp();
    my $line = $_;
        
    my ($target_id, $est_counts) = split(/\t/, $line) ;

    if( $line_number >= 0)   # \s == wide space; \s+ == any space; [0-9]+ == any integer, any times; ([0-9]+) == put in special brackets so pts in into special variable $1. if (()) then $2.
      {
        $hash_ENST{$target_id} +=  $est_counts;
      } 
      $line_number++;
 }


close(FILE2);


#open(OUT, ">$ARGV[2]") || die "Can't open file for writing: $!\n"; # > ==write to the file; >> == append to file
#print "target_id\test_counts\n";

# $i is a key here
# $hash_ENST{$i} is value
my $cnt = 0;
foreach my $i(sort keys %hash_ENST)
  {
    
    printf "$i\t$hash_ENST{$i}\n" if($cnt >= 0);
    #printf OUT "$i\t$hash_ENST{$i}\n";

    $cnt++;
  }

#close(OUT);
#print "job done\n";





