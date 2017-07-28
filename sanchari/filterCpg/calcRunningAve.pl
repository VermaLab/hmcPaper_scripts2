#!/usr/bin/perl

use strict;
use warnings;

my $N = 1;
my $curAve = 0;
while (<STDIN>){
	my @fields = split '#';
	my ($a, $b) = split ' ', $fields[1];
	#TODO:  Can we not look at less than 100 percent methylation. 
	#TODO:  Like maybe 50 and above?
	$curAve = $curAve + (($a/$b) - $curAve)/$N;
	$N++;
	if ($N % 100000 == 0){
		#print "$a\t$b\t$curAve\t$N\n";
		print STDERR "$curAve\t$N\n";
	}
}

print "$curAve\t$N\n";
