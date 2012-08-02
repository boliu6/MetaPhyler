#!/usr/bin/perl

#############################################
#
# Program: Merge multiple taxprof files in to one table
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
if (scalar @ARGV == 0) {
    Usage();
}
#----------------------------------------#


my %tax = ();
my %sams = ();
foreach my $file (@ARGV) {
    $file =~ /^(\S+?)\./;
    my $sam = $1;
    $sams{$sam} = 1;
    
    open(FH, "$file");
    foreach my $line (<FH>) {
	if ($line =~ /^Name/) {
	    next;
	}
	my ($id, $pct) = split("\t", $line);
	$tax{$id}{$sam} = $pct;
    }
    close FH;
}

print "\t", join("\t", sort keys %sams), "\n";
foreach my $tax (sort keys %tax) {
    print "$tax\t";
    foreach my $sam (sort keys %sams) {
	if (exists $tax{$tax}{$sam}) {
	    print "$tax{$tax}{$sam}\t";
	}
	else {
	    print "0.00\t";
	}
    }
    print "\n";
}

exit;


sub Usage {
    die("
Usage:
       perl report.pl <taxprof 1> <taxprof 2> ...

       Merge multiple taxonomy profile files into one table.

Options:
       <taxprof>      Taxonomy profile output from program taxprof

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
