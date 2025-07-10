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

        Compare case with baseline mean and zscore.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 case1.pfe baseline.OLOZscore prefix

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $f = shift;
my $bf = shift;
my $prefix = shift;

my %bl;
open IN,"$bf" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = shift @t;
	$bl{$id} = \@t;
}
close IN;



my %gh;

&calScore($f,\%gh);

#print Dumper(\%gh);
# unique to the case
open CASE,">","./$prefix.spcific.feature.txt" || die $!;
foreach my $k (keys %gh){
	unless(exists $bl{$k}){
		my @k = split /#/,$k;
		my $score = $gh{$k}{"s"};
		my $vscore = $gh{$k}{"sv"};
		
		if($score > 0.2 && abs($vscore) < 0.2){
			print CASE join "\t",@k;
			print CASE "\t$score";
			print CASE "\t$vscore\n";
		}
	}
}
close CASE;

# elevated regions
open OUT,">","./$prefix.elvated.lst" || die $!;
foreach my $k (keys %gh){
	if(exists $bl{$k}){
		my @k = split /#/,$k;
		my $arr = $bl{$k};
		my $score = $gh{$k}{"s"};
		my $vscore = $gh{$k}{"sv"};
		my @zarr;
		for(my $i = 0; $i < scalar(@$arr); $i=$i+2){
			my $j = $i + 1;
			my $m = $arr->[$i];
			my $sd = $arr->[$j];
			my $z = ($score-$m)/$sd;
			push @zarr,$z;
		}

		my $lcount = 0;
		my $rcount = 0;
		my $mcount = 0;
		foreach my $e (@zarr){
			if($e >= 1.96){
				$rcount++;
			}elsif($e <= -1.96){
				$lcount++;
			}else{
				$mcount++;
			}
		}

		# setting criteria for output
		if($rcount/scalar(@zarr) >= 0.9 && $score > 0.2){
			print OUT join "\t",@k;
			print OUT "\t";
			print OUT join "\t",@zarr;
			print OUT "\n";
		}
	}
}
close OUT;
#my $i = 0;
#foreach my $gf (@g){
#	my %gh;
#	&calScore($gf,\%gh);
#	
#	foreach my $k (keys %gh){
#		if(exists $gh{$k}){
#			my $score = $gh{$k}{"s"};
#			$sindex{$k}->[$i] = $score;
#		}
#	}
#	$i++;
#	%gh = ();
#}
#
#$i = 0;
#foreach my $gf (@g){
#	my %gh;
#	&calScore($gf,\%gh);
#	
#	foreach my $k (keys %gh){
#		if(exists $gh{$k}){
#			my $score = $gh{$k}{"sv"};
#			$svindex{$k}->[$i] = $score;
#		}
#	}
#	$i++;
#	%gh = ();
#}
#
##print Dumper(\%sindex);
#
## score
#foreach my $k (keys %sindex){
#	my $str_tmp = join "",@{$sindex{$k}};
#	if($str_tmp =~ /NA/){
#		next;
#	}else{
#		my $arr = $sindex{$k};
#		my @oarr;
#		for(my $i = 0; $i < scalar(@$arr); $i++){
#			my @tarr;
#			for(my $j = 0; $j < scalar(@$arr); $j++){
#				if($i != $j){
#					push @tarr,$arr->[$j];
#				}
#			}
#			my @r;
#			calOneLeaveOutZscore(\@tarr,\@r);
#			push @oarr,@r;
#		}
#		
#		print "$k\t";
#		print join "\t",@oarr;
#		#print "\t";
#		#print join "\t",@$arr;
#		print "\n";
#	}
#}
#
#sub calOneLeaveOutZscore {
#	my $larr = shift;
#	my $zarr = shift;
#	my $lstd = std_dev(@$larr);
#	my $lmean = average(@$larr);
#	push @$zarr,$lmean;
#	push @$zarr,$lstd;
#	#my $sample_count = scalar(@$larr);
#	#my $lmean = sum(@$larr)/$sample_count;
#	
#}
#
#sub average {
#        my (@values) = @_;
#
#        my $count = scalar @values;
#        my $total = 0; 
#        $total += $_ for @values; 
#
#        return $count ? $total / $count : 0;
#}
#
#sub std_dev {
#        my ($average, @values) = @_;
#
#        my $count = scalar @values;
#        my $std_dev_sum = 0;
#        $std_dev_sum += ($_ - $average) ** 2 for @values;
#
#        return $count ? sqrt($std_dev_sum / $count) : 0;
#}

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
		}elsif($nan_count <= 2){
			# replace NaN
			my @clean;
			foreach my $pfe (@t){
				if($pfe ne "nan"){
					push @clean,$pfe;
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
			my $vscore = 1 - abs(($score+0.01)/(2*($mid_value-$boundary_value)+0.01));
			#if($score < 0){
			#	$score = 0;
			#}
			#if($vscore eq ""){
			#	print STDERR "warning\n";
			#}
			$h->{$id}{"s"} = $score;
			$h->{$id}{"sv"} = $vscore;
			#print STDERR "$chr\t$s\t$e\t$n1\t$n2\t$strand\t$score\t$vscore\n";
		}
	}
	close IN;
}

print STDERR "analysis done\n";
print STDERR `date`;
