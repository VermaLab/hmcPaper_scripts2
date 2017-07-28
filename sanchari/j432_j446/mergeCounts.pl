#!/usr/bin/perl

use strict;
use warnings;

#load in the command line arguments
my ($bed1, $bed2, $chr) = @ARGV;

#print "file1:  $bed1\n";
#print "file2:  $bed2\n";
#print "chr:  $chr\n";



#open file handles to each file
open my $in1, '<', $bed1 or die;
open my $in2, '<', $bed2 or die;


#read the first line in each file and init the fields
my ($chr1, $pos1, $nMeth1, $nTot1) = split ' ', <$in1>;
my ($chr2, $pos2, $nMeth2, $nTot2) = split ' ', <$in2>;

#print "$chr1, $pos1, $nMeth1, $nTot1\n";
#print "$chr2, $pos2, $nMeth2, $nTot2\n";

#keep going while there's still stuff
#being read by both files
while (!eof $in1 and !eof $in2){
	#same position
	#add counts from both lines
	#move both file handles up one
	if ($pos1 == $pos2){
		print "$chr1\t$pos1\t".($nMeth1+$nMeth2)."\t".($nTot1+$nTot2)."\n";
		($chr1, $pos1, $nMeth1, $nTot1) = split ' ', <$in1>;
		($chr2, $pos2, $nMeth2, $nTot2) = split ' ', <$in2>;
	}
	elsif ($pos1 > $pos2){
		#print 2nd file, move it's handle up one 	
		print "$chr2\t$pos2\t$nMeth2\t$nTot2\n";
		($chr2, $pos2, $nMeth2, $nTot2) = split ' ', <$in2>;
	}
	elsif ($pos2 > $pos1){
		#print 1st file, move it's handle up one 	
		print "$chr1\t$pos1\t$nMeth1\t$nTot1\n";
		($chr1, $pos1, $nMeth1, $nTot1) = split ' ', <$in1>;
	}
	#print "1:".<$in1>;
	#print "2:".<$in2>;
}

#one of the files finished, now write out the rest
if (eof $in1){
	while(<$in2>){
		print "$chr2\t$pos2\t$nMeth2\t$nTot2\n";
		($chr2, $pos2, $nMeth2, $nTot2) = split ' ', <$in2>;
	}
}
else{
	while(<$in1>){
		print "$chr1\t$pos1\t$nMeth1\t$nTot1\n";
		($chr1, $pos1, $nMeth1, $nTot1) = split ' ', <$in1>;
	}
}
