#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        Design for get tss from a gene loci bed file.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 XX.bed

USAGE
print "$usage";
exit(1);
};

my $bed=shift;

open IN,"$bed" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	
	my $strand = $t[-1];
	my ($s,$e);
	if($strand eq "+"){
		$s = $t[1];
		$e = $s + 1;
	}elsif($strand eq "-"){
		$s = $t[2]-1;
		$e = $t[2]
	}else{
		$s =  int(($t[1] + $t[2])/2);
		$e = $s + 1;
	}
	print "$t[0]\t$s\t$e\t$t[3]\t$t[4]\t$t[5]\n";
}
close IN;

