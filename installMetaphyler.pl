#!/usr/bin/perl

#############################################
#
# Program: Install Metaphyler programs.
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Sat Jul 28 14:37:59 EDT 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


# format blast database
my $cmd = "formatdb -p F -i $Bin/markers/markers.dna";
print "$cmd\n";
system($cmd);

$cmd = "formatdb -p T -i $Bin/markers/markers.protein";
print "$cmd\n";
system($cmd);

my $gcc = "g++ -Wall -W -O2";
my @programs = ("simuReads", "metaphylerClassify", "taxprof", "combine", "scores");
foreach my $program (@programs) {
    $cmd = "$gcc -o $Bin/bin/$program $Bin/src/$program.cpp";
    print "$cmd\n";
    system($cmd);
}

exit;
