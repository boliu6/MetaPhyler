#!/usr/bin/perl

#############################################
#
# Program: Classify sequences.
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Fri Mar  9 23:24:30 EST 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);

#----------------------------------------#
# read command line options
#----------------------------------------#
my $query = "";
my $blast = "";
my $prefix = "";
my $nump = 0;
if (scalar @ARGV == 4) {
    ($query, $blast, $prefix, $nump) = @ARGV;
    if ($blast ne "blastn" && $blast ne "blastx") { Usage();}
} else {
    Usage();
}
#----------------------------------------#


my $ref = "$Bin/markers/markers.dna";
my $param = "-W15";
if ($blast eq "blastx") {
    $param = "";
    $ref = "$Bin/markers/markers.protein";
}
# run blast
my $cmd = "blastall -p $blast $param -a$nump -e0.01 -m8 -b1 -i $query -d $ref > $prefix.$blast";
print "$cmd\n";
system("$cmd");

# classification
$cmd = "$Bin/metaphylerClassify $Bin/markers/markers.$blast.classifier $Bin/markers/markers.taxonomy $prefix.$blast > $prefix.classification";
print "$cmd\n";
system("$cmd");

$cmd = "$Bin/taxprof 0.9 $prefix.classification $prefix $Bin/markers/tid2name.tab";
print "$cmd\n";
system("$cmd");

exit;


sub Usage {
    die("
Usage:
       perl runMetaphyler.pl <query> <blast> <prefix> <# threads>

Options:
       <query>        Query sequences in FASTA format to be classified.
       <blast>        blastn or blastx.
                      You can try both modes, and combine the classifications. 
                      blastn is recommended for short reads (100bp).
       <prefix>       Output prefix.
       <# threads>    Number of threads to run BLAST.

Output:
       prefix.blast[n/x]
                      Raw blast output.

       prefix.classification
                      Classification results.

       prefix.<genus|family|order|class|phylum>.taxprof
                      Taxonomy profiles at each level.

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
