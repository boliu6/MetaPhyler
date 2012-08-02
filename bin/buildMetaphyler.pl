#!/usr/bin/perl

#############################################
#
# Program: Perform simulation, build metaphyler
#          classification models.
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
my $norm = "";
my $seqtype = "";
my $lens = "";
my $step = 30;
my $qfile = "";
my $rfile = "";
my $taxfile = "";
my $pre = "";
my $nump = 0;
my $blast = "";
if (scalar @ARGV == 8) {

    if    ($ARGV[0] eq "norm") { $norm = "true";}
    else                       { $norm = "false";}

    ($qfile, $rfile, $lens, $taxfile, $blast, $pre, $nump) = @ARGV[1...7];
} else {
    Usage();
}
#----------------------------------------#

# format blast database and run blast
my $p = "F";
my $param = "-FF -W15 -a$nump";
if ($blast eq "blastx" || $blast eq "blastp") {
    $p = "T";
    $param = "";
}
my $cmd = "formatdb -p $p -i $rfile";
print "$cmd\n";
system("$cmd");

my @lens = split(',', $lens);
my $outfile = "$pre.$blast.classifier";
if (-e $outfile) {
    system("rm $outfile");
}

foreach my $len (@lens) {

    my $prefix = "$pre.$len";
# simulate reads
    my $cmd = "$Bin/simuReads $len $step $qfile > $prefix.fasta";
    print "$cmd\n";
    system("$cmd");
    
$cmd = "blastall -p $blast $param -e1e-3 -m8 -b1000 -v1000 -i $prefix.fasta -d $rfile > $prefix.$blast";
    print "$cmd\n";
    system("$cmd");
    
# train model
    $cmd = "$Bin/metaphylerTrain norm $taxfile $rfile $prefix.$blast $len $blast >> $pre.$blast.classifier";
    print "$cmd\n";
    system("$cmd");
    
}
exit;


sub Usage {
    die("
Usage:
       perl buildMetaphyler.pl <norm|unnorm> <fasta 1> <fasta 2> <lengths> <taxonomy> <blast> <prefix> <# threads>

Options:
       <norm|unnorm>  Perform normalization (true) or not (false).
                      If the taxonomy is well defined, then normalization is recommended.
       <lengths>      Lengths of simulated reads, based on which classifiers are built.
                      To build classifiers for 200bp, 300bp and 400bp, use 200,300,400.
       <fasta 1>      From which reads are simulated.
       <fasta 2>      To which simulated reads are mapped.
                      For blastn, <fasta 1> and <fasta 2> are the same file.
       <taxonomy>     Taxonomy labels of sequences in <fasta 2>.
       <prefix>       Output prefix.
       <# threads>    Number of threads to run BLAST.

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
