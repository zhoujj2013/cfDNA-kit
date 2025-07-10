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

        Analysis cfDNA fragment length diversity (measured by PFE) in given regions in a specific sample.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xx.bam xx.region.tss.bed xx.chrom.size.txt prefix

USAGE
print "$usage";
exit(1);
};

print `date`;
print STDERR "analysis start\n";

my $bam = shift;
my $tss_tmp = shift;
my $gsize = shift;
my $prefix = shift;

# filter the region will exceed the boundary of the genome
# add this part when running in rabbit genome, because the fragment of the genome is short.
#
my %gsize;
open IN,"$gsize" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$gsize{$t[0]} = $t[1];
}
close IN;

open BED,">","$prefix.filtered.tss.bed" || die $!;
open IN,"$tss_tmp" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	if(!exists $gsize{$t[0]}){ # remove the gene not in ref. genome
		next;
	}elsif($t[1] - 5000 <= 0){ # remove the gene in left boundary
		next;
	}elsif($t[2] + 5000 >= $gsize{$t[0]}){ # remove the gene in the right boundary
		next;
	}else{  # output the remaining genes.
		print BED "$_\n";
	}
}
close IN;
close BED;

my $tss = "$prefix.filtered.tss.bed";

# create cell files
open OUT,">","run.$$.sh" || die $!;

print OUT "python3 $Bin/cal_pfe.py $bam $tss -5000 -4000 > $prefix.m5000.m4000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss -4000 -3000 > $prefix.m4000.m3000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss -3000 -2000 > $prefix.m3000.m2000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss -2000 -1000 > $prefix.m2000.m1000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss -1000 0 > $prefix.m1000.0.pfe\n";

print OUT "python3 $Bin/cal_pfe.py $bam $tss 0 1000 > $prefix.0.1000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss 1000 2000 > $prefix.1000.2000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss 2000 3000 > $prefix.2000.3000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss 3000 4000 > $prefix.3000.4000.pfe\n";
print OUT "python3 $Bin/cal_pfe.py $bam $tss 4000 5000 > $prefix.4000.5000.pfe\n";

close OUT;
`perl $Bin/multi-process.pl -cpu 10 run.$$.sh`;

# merge the files
my @f = ("$prefix.m5000.m4000.pfe", "$prefix.m4000.m3000.pfe", "$prefix.m3000.m2000.pfe", "$prefix.m2000.m1000.pfe","$prefix.m1000.0.pfe","$prefix.0.1000.pfe","$prefix.1000.2000.pfe","$prefix.2000.3000.pfe","$prefix.3000.4000.pfe","$prefix.4000.5000.pfe");

my %loc2pfe;
foreach my $ff (@f){
	open IN,"$ff" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		pop(@t);
		my $id = join "#",@t;
		my @arr = ("nan","nan","nan","nan","nan","nan","nan","nan","nan","nan");
		$loc2pfe{$id} = \@arr;
	}
	close IN;
}

for(my $i = 0; $i < scalar(@f); $i++){
	open IN,"$f[$i]" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $pfe = pop(@t);
		my $id = join "#",@t;
		$loc2pfe{$id}->[$i] = $pfe; 
	}
	close IN;
}

open OUT,">","$prefix.merged.pfe" || die $!;
foreach my $id (keys %loc2pfe){
	my $outarr = $loc2pfe{$id};
	my @id = split /#/,$id;
	print OUT join "\t",@id;
	print OUT "\t";
	print OUT join "\t",@$outarr;
	print OUT "\n";
}
close OUT;

`Rscript  ~/github/jjUtil/compareDisRNAseq/MultiBoxplotNoOutlineSize.r $prefix 0.5 6 6 $prefix.m5000.m4000.pfe:7:gray $prefix.m4000.m3000.pfe:7:gray $prefix.m3000.m2000.pfe:7:gray $prefix.m2000.m1000.pfe:7:gray $prefix.m1000.0.pfe:7:gray $prefix.0.1000.pfe:7:gray $prefix.1000.2000.pfe:7:gray $prefix.2000.3000.pfe:7:gray $prefix.3000.4000.pfe:7:gray $prefix.4000.5000.pfe:7:gray`;

print STDERR "Rscript  ~/github/jjUtil/compareDisRNAseq/MultiBoxplotNoOutlineSize.r $prefix 0.5 6 6 $prefix.m5000.m4000.pfe:7:gray $prefix.m4000.m3000.pfe:7:gray $prefix.m3000.m2000.pfe:7:gray $prefix.m2000.m1000.pfe:7:gray $prefix.m1000.0.pfe:7:gray $prefix.0.1000.pfe:7:gray $prefix.1000.2000.pfe:7:gray $prefix.2000.3000.pfe:7:gray $prefix.3000.4000.pfe:7:gray $prefix.4000.5000.pfe:7:gray\n\n";

print STDERR "analysis done\n";
`rm run.$$.sh `;
my $fstr = join " ",@f;
`rm $fstr`;
print `date`;
