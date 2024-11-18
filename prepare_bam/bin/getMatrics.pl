#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
#use lib "/home/zhoujj/my_lib/pm";
#use bioinfo;

&usage if @ARGV<1;

#open IN,"" ||die "Can't open the file:$\n";
#open OUT,"" ||die "Can't open the file:$\n";

sub usage {
        my $usage = << "USAGE";

        Design for get statistical matrics for ATAC-seq analysis.
        Author: zhoujj2013\@gmail.com
        Usage: $0 sample_prefix1 sample_prefix2 sample_prefix3 ...

USAGE
print "$usage";
exit(1);
};

# for triming, all are fragment
print "Trim low quanlity and adaptor\n";
print "\tFiltered_flagment\tTotal_flagment\tPassed_flagment\tPassRate\n";
my @readsfiltered;
my @readstotal;
my @readspass;
my @readsper;
foreach my $s (@ARGV){
	my @f = glob("$s/00datafilter/*adapterRemoval.setting");
	my $f = $f[0];
	&checkExists($f);
	my $total = 0;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		next unless(/Total number of read pairs:/);
		$total = $1 if(/Total number of read pairs: (\S+)/);
	}
	close IN;
	
	my $pass = 0;
	my @f1 = glob("$s/00datafilter//$s.adapterRemoval.setting"); 
	my $f1 = $f1[0];
	&checkExists($f1);
	open IN,"$f1" || die $!;
	while(<IN>){
		chomp;
		next unless(/Number of retained reads:/);
		$pass = $1 if(/Number of retained reads:\s+(\d+)/);
		$pass = int($pass/2);
	}
	close IN;
	
	my $filtered = $total - $pass;
	my $per = sprintf("%.2f",($pass)*100/$total);
	
	push @readsfiltered,$filtered;
	push @readstotal,$total;
	push @readspass,$pass;
	push @readsper,$per;	
}
ArrAsCol(\@ARGV,\@readsfiltered,\@readstotal,\@readspass,\@readsper);
print "\n";

#for mapping rate
print "Alignment Summary\n";
print "\tTotal_flagment\tAligned_flagment\tAlignmentRate(%)\tproperlyPaired\tproperlyPairedRate\tdup\tdupRate\tUnmappedFlag\n";
my @AlignInput;
my @aligned;
my @AlignRate;
my @AlignedProperly;
my @AlignedProperlyRate;
my @DupCount;
my @DupRate;
my @UnmappedFlag;
my $i = 0;
foreach my $s (@ARGV){
	my $f = "$s/01alignment/$s.BackgroudAlign.sorted.bam.flagstat";
	&checkExists($f);
	my $readCleanTotal;
	my $aligned_count;
	my $aligned_rate;
	my $aligned_properly;
	my $aligned_properly_rate;

	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		#next unless(/concordant pair alignment rate/);
		$readCleanTotal = $1 if(/(\d+)\s\+\s\d+\sin total \(/); # 977065807 + 0 in total (QC-passed reads + QC-failed reads)
		$aligned_count = $1 if(/(\d+)\s\+\s\d+\smapped \(/); # 975017234 + 0 mapped (99.79% : N/A)

		$aligned_properly = $1 if(/(\d+)\s\+\s\d+\sproperly paired \(/); #834624816 + 0 properly paired (85.91% : N/A)
		#$uniq_aligned = $1 if(/Uniquely mapped reads number\s+\|\s+(\d+)/);
		#$multi_aligned = $1 if(/Number of reads mapped to multiple loci\s+\|\s+(\d+)/);
		#$rate = $1 if(/(\S+)% concordant pair alignment rate/);
		#$aligned = int($readCleanTotal * $rate * 0.01);
	}
	close IN;
	$readCleanTotal = int($readCleanTotal/2);
	$aligned_rate = sprintf "%.2f",(($aligned_count/2)/$readCleanTotal)*100;
	$aligned_properly_rate = sprintf "%.2f",(($aligned_properly/2)/$readCleanTotal)*100;
	push @AlignInput,$readCleanTotal;
	push @AlignRate,$aligned_rate;
	push @aligned,int($aligned_count/2);
	push @AlignedProperly,int($aligned_properly/2);
	push @AlignedProperlyRate,$aligned_properly_rate;


	my $f1 = "$s/01alignment/$s.rmdup.log";
	&checkExists($f1);
	my $dup_rate;
	my $dup_count;
	open IN,"$f1" || die $!;
	while(<IN>){
		chomp;
		next unless(/^Unknown/);
		my @t = split /\t/;
		$dup_rate = $t[8];
	}
	close IN;
	$dup_count = int($aligned_count*$dup_rate);
	push @DupCount,int($dup_count/2);
	push @DupRate,$dup_rate*100;

	my $f2 = "$s/01alignment/$s.unmappedR1.fastq";
	&checkExists($f1);
	my $unmapped_flag = `wc -l $f2 | awk '{print \$1}'`;
	chomp($unmapped_flag);
	push @UnmappedFlag,$unmapped_flag;
		
	$i++;
}
$i=0;
ArrAsCol(\@ARGV,\@AlignInput,\@aligned,\@AlignRate,\@AlignedProperly,\@AlignedProperlyRate,\@DupCount,\@DupRate,\@UnmappedFlag);
print "\n";

#for combine insert pdfs
my @pdfs;

foreach my $s (@ARGV){
	my $pdf = "$s/01alignment/$s.insert.hist.pdf";
	push @pdfs,$pdf;
}
my $pdfs_str = join " ",@pdfs;

`convert -quality 100 -density 100 $pdfs_str merged.pdf`;

print "Please check merged.pdf for insert information\n";


sub ArrAsCol{
	my $arr = \@_;
	my $len = scalar(@{$arr->[0]});
	for(my $i = 0; $i < $len; $i++){
		my @l;
		foreach my $a (@$arr){
			push @l,$a->[$i];
		}
		print join "\t",@l;
		print "\n";
	}
}

sub checkExists{
	my $ffname = shift @_;
	unless(-e $ffname){
		print STDERR "$ffname not exists!\n";
		exit(1);
	}
}
