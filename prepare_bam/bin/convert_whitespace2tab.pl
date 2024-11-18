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

        This script create makefile for xxx.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

my $f=shift;
open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my $l = $_;
	$l =~ s/^(\s+)//g;
	my @t = split /\s+/,$l;
	my $count = shift @t;
	my $id = join " ",@t;
	#print scalar(@t),"\n";
	print "$count\t$id\n";
}
close IN;
