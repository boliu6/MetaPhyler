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
my $ref = "";
my $taxfile = "";
my $model = "";
my $blast = "";
my $prefix = "";
my $nump = 0;
if (scalar @ARGV == 7) {
    ($query, $ref, $taxfile, $model, $blast, $prefix, $nump) = @ARGV;
    if ($blast ne "blastn" && $blast ne "blastp" && $blast ne "blastx") { Usage();}
} else {
    Usage();
}
#----------------------------------------#


# format blast database and run blast
my $param = "-FF -a$nump";
if ($blast eq "blastn") {
    $param = "-FF -W20 -a$nump";
}

my $cmd = "blastall -p $blast $param -e1e-3 -m8 -b1 -v1 -i $query -d $ref > $prefix.$blast";
print "$cmd\n";
system("$cmd");

# classification
$cmd = "$Bin/metaphylerClassify $model $taxfile $prefix.$blast > $prefix.$blast.classification";
print "$cmd\n";
system("$cmd");


exit;


sub Usage {
    die("
Usage:
       perl runMetaphyler.pl <query> <reference> <taxonomy> <classifiers> <blast> <prefix> <# threads>

Options:
       <query>        Input sequences to be classified.
       <reference>    Known sequences with taxonomy lables.
       <taxonomy>     Taxonomy labels of the reference.
       <classifiers>  Models built from buidMetaphyler.pl.
                      If have multiple models, separate them with comma (e.g., fileA,fileB).
                      Example: fileA is built from 100bp, and fileB is built from 200bp.
                      You can also put them in to one file.
       <blast>        blastn, blastp or blastx.
       <prefix>       Output prefix.
       <# threads>    Number of threads to run BLAST.

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
