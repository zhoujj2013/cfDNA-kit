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
        Usage: $0 samples.lst

USAGE
print "$usage";
exit(1);
};

my $conf=shift;

my @pre_index;
my %total;
my %index;
open IN,"$conf" || die $!;
while(<IN>){
	chomp;
	next if(/^#/);
	my @t = split /\t/;
	my $id = $t[0];
	push @pre_index,$id;
	open STAT,"out/$id/01alignment/$id.BackgroudAlign.sorted.bam.flagstat" || die $!;
	while(<STAT>){
		chomp;
		next unless(/total/);
		my $c = $1 if(/(\d+)/);
		$total{$id} = $c;
	}
	close STAT;
	
	open STAT,"out/$id/02kraken/$id.spe.count.stat.rf" || die $!;	
	while(<STAT>){
		chomp;
		my @cc = split /\t/;
		$index{$cc[1]} = 1;
	}
	close STAT;
}
close IN;


my %out;
open IN,"$conf" || die $!;
while(<IN>){
	chomp;
	next if(/^#/);
	my @t = split /\t/;
	my $id = $t[0];

	my %s;
	open STAT,"out/$id/02kraken/$id.spe.count.stat.rf" || die $!;
	while(<STAT>){
		chomp;
		my @cc = split /\t/;
		$s{$cc[1]} = $cc[0];
	}
	close STAT;
	
	my $total_c = $total{$id};
	
	foreach my $k (sort keys %index){
		if(exists $s{$k}){
			# read count per 10M
			my $abd = $s{$k}/($total_c/10000000);
			push @{$out{$id}},[$abd,$s{$k},$total_c];
		}else{
			my $abd = 0/($total_c/10000000);
			push @{$out{$id}},[$abd,0,$total_c];
		}
	}
}
close IN;

my @index_sorted = sort keys %index;

print "geneid\t";
print join "\t", @pre_index; #sort keys %out;
print "\t";
print join "\t", @pre_index;
print "\t";
print join "\t", @pre_index;
print "\n";

#print "group";
#foreach my $sample_id (sort keys %out){
#	my @ss = split /_/,$sample_id;
#	print "\t$ss[2]";
#}
#print "\n";

#print Dumper(\%out);

for(my $i = 0; $i < scalar(@index_sorted); $i++){
	my @s2c;
	my @s2p;
	my @s2t;
	#my $flag = 0;
	my $nonzeroCount = 0;
	foreach my $sid (@pre_index){ #sort keys %out){
		my $count = $out{$sid};
		push @s2c,$count->[$i]->[1];
		push @s2p,$count->[$i]->[0];
		push @s2t,$count->[$i]->[2];
		if($count->[$i]->[1] >= 1){
			$nonzeroCount++;
		}
	}
	#if($nonzeroCount > 1){
		print "$index_sorted[$i]\t";
		print join "\t",@s2c;
		print "\t";
		print join "\t",@s2p;
		print "\t";
		print join "\t",@s2t;
		print "\n";
	#}
}

#print Dumper(\%out);
