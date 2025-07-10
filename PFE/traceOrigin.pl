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

        Calculate the correlation of detal pfe and gene expression of different tissue and cell types.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 tissue.expr.mat pfe.peaks.lst > xx.cor
        column: tissue
        row: sample_name

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $exprf = shift;
my $peakf = shift;

my %e;
my %g;
open IN,"$exprf" || die $!;

my $header = <IN>;
chomp($header);
my @h = split /\t/,$header;

while(<IN>){
	chomp;
	my @t = split /\t/;
	my $geneid = $t[0];
	for(my $i = 1; $i < scalar(@t); $i++){
		$e{$h[$i]}{$geneid} = $t[$i];
		$g{$geneid} = 1;
	}
}
close IN;

#print Dumper(\%e);

open IN,"$peakf" || die $!;

$h[0] = "Sample_id";
print join "\t",@h;
print "\n";

while(<IN>){
	chomp;
	my @t = split /\t/;
	my $prefix = shift @t;
	my $f = shift @t;
	my %case;

	open PEAK,"$f" || die $!;
	while(<PEAK>){
		chomp;
		my @t = split /\t/;
		my $delta_pfe = $t[6];
		my $geneid = $t[3];
		if(exists $g{$geneid}){
			$case{$geneid} = $delta_pfe;
		}
	}
	close PEAK;
	
	my @rho_arr;
	for(my $j = 1; $j < scalar(@h); $j++){
		# all zetero
		my @expr_arr;
		open OUT,">","tmp-$j.txt" || die $!;
		foreach my $geneid (sort keys %case){
			print OUT "$case{$geneid}\t$e{$h[$j]}{$geneid}\n";
			my $ee = $e{$h[$j]}{$geneid};
			push @expr_arr,$ee;
		}
		close OUT;
		my $sum_expr = sum(@expr_arr);
		if($sum_expr == 0){
			push @rho_arr,"0";
		}else{
			my $result = `Rscript $Bin/cal_cor.r ./tmp-$j.txt`;
			my @result = split /\s+/,$result;
			my $rho = sprintf ("%.4f", $result[1]);
			my $p = $result[2];
			push @rho_arr,$rho;
		}
		`rm ./tmp-$j.txt`;
	}
	print "$prefix\t";
	print join "\t",@rho_arr;
	print "\n";
}
close IN;


print STDERR "analysis done\n";
print STDERR `date`;
