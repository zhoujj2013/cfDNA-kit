#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
use List::Util qw(sum);
use Statistics::Basic qw(:all);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        Convert correlation to ranking.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xx.cor > xx.cor.rank
        column: tissue
        row: sample_name

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $corf = shift;
open IN,"$corf" || die $!;

my $header = <IN>;
chomp($header);
my @header = split /\t/,$header;
print join "\t",@header;
print "\n";

while(<IN>){
	chomp;
	my @t = split /\t/;
	#my $id = shift @t;
	my %c;
	for(my $i = 1; $i < scalar(@t); $i++){
		$c{$header[$i]} = $t[$i];
	}

	my @crank = sort{$c{$a} <=> $c{$b}} keys %c;
	
	my %crank;
	for(my $j = 0; $j < scalar(@crank); $j++){
		my $tname = $crank[$j];
		$crank{$tname} = $j + 1;
	}

	for(my $i = 1; $i < scalar(@t); $i++){
		my $tname = $header[$i];
		$t[$i] = $crank{$tname};
	}
	print join "\t",@t;
	print "\n";
}
close IN;

print STDERR "analysis done\n";
print STDERR `date`;
