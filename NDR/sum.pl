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

        This script create makefile for LncFunNet analysis.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

my $f = shift;

open IN,"$f" || die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $name = shift @t;
	my $sum = 0;
	foreach my $i (@t){
		$sum = $sum + $i;
	}
	print "$name\t$sum\n";
}
close IN;
