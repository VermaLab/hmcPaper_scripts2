#!/usr/bin/perl
#print the chromosomal locations that have at least k reads
#read a conversionCount file from stdin
#write locations to stdout
#usage:
#   zcat HPNE_BS_3.txt.gz | ./filterSitesByCount.pl 4 > common4_HPNE_BS_3.txt
use strict;
use warnings;

my $k = 0;
if (scalar @ARGV > 0){
	$k = $ARGV[0];
	#print "k:  $k\n";
}else{
	die "must supply read count threshold argument\n";
}

while(<STDIN>){
	#my $line = $_;
	my @vals = split;
	if ($vals[3] >= $k){
		#print $line;
		print "$vals[0]\t$vals[1]\n";
	}
}
