#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
use List::Util qw(sum);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        Prerank genes by delta pfe value (score value).
        Author: zhoujj2013\@gmail.com 
        Usage: $0 case1.pfe

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $f = shift;

my %gh;

&calScore($f,\%gh);

#print Dumper(\%gh);
# unique to the case
foreach my $k (keys %gh){
	my @k = split /#/,$k;
	my $score = $gh{$k}{"s"};

	print join "\t",@k;
	print "\t$score\n";
	
}


# read in the files
sub calScore {
	my $f = shift;
	my $h = shift;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $chr = shift @t;
		my $s = shift @t;
		my $e = shift @t;
		my $n1 = shift @t;
		my $n2 = shift @t;
		my $strand = shift @t;
		
		my $id = "$chr#$s#$e#$n1#$n2#$strand";
		
		# deal with NaN
		my $nan_count = 0;
		foreach my $pfe (@t){
			if($pfe eq "nan"){
				$nan_count++;
			}
		}
		if($nan_count > 2){
			next;
		}
		# replace NaN
		my @clean;
		if($nan_count <= 2){
			foreach my $pfe (@t){
				if($pfe ne "nan"){
					push @clean,$pfe;
				}
			}
		}
		my $mean = sum(@clean)/scalar(@clean);
		for(my $i = 0; $i < scalar(@t); $i++){
			if($t[$i] eq "nan"){
				$t[$i] = $mean;
			}
		}
	
		# calculate peak point
		my $point_value = ($t[4]+$t[5])/2;
		my $boundary_value = ($t[0]+$t[1]+$t[8]+$t[9])/4;
		my $mid_value = ($t[2]+$t[3]+$t[6]+$t[7])/4;
	
		my $score = $point_value - $boundary_value;
		my $vscore = 1 - abs($score/(2*($mid_value-$boundary_value)));
		#if($score < 0){
		#	$score = 0;
		#}
		$h->{$id}{"s"} = $score;
		$h->{$id}{"sv"} = abs($vscore);
		#print "$chr\t$s\t$e\t$n1\t$n2\t$strand\t$score\t$vscore\n";
	}
	close IN;
}

print STDERR "analysis done\n";
print STDERR `date`;
