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

        Build baseline for a normal group. Prepare for comparision with a disease sample.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 s1.pfe,s2.pfe...

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $g1_str = shift;

my @g1 = split /,/,$g1_str;

# build the index
my @g = @g1;

#####################################
my @na;
foreach my $g (@g){
	push @na,"NA";
}
#####################################

my %sindex;
foreach my $gf (@g){
	my %gh;
	&calScore($gf,\%gh);
	foreach my $k (keys %gh){
		$sindex{$k} = [];
		#if(exists $sindex{$k}){
		#	push @{$sindex{$k}},"NA";
		#}else{
		#	$sindex{$k} = [];
		#	push @{$sindex{$k}},"NA";
		#}
	}
}

foreach my $k (keys %sindex){
	my @na_tmp = @na;
	$sindex{$k} = \@na_tmp;
}

my %svindex;
foreach my $gf (@g){
        my %gh;
        &calScore($gf,\%gh);
        foreach my $k (keys %gh){
		$svindex{$k} = [];
                #if(exists $svindex{$k}){
                #        push @{$svindex{$k}},"NA";
                #}else{
                #        $svindex{$k} = [];
                #        push @{$svindex{$k}},"NA";
                #}
        }
}

foreach my $k (keys %svindex){
        my @na_tmp = @na;
        $svindex{$k} = \@na_tmp;
}


my $i = 0;
foreach my $gf (@g){
	my %gh;
	&calScore($gf,\%gh);
	
	foreach my $k (keys %gh){
		if(exists $gh{$k}){
			my $score = $gh{$k}{"s"};
			$sindex{$k}->[$i] = $score;
		}
	}
	$i++;
	%gh = ();
}

$i = 0;
foreach my $gf (@g){
	my %gh;
	&calScore($gf,\%gh);
	
	foreach my $k (keys %gh){
		if(exists $gh{$k}){
			my $score = $gh{$k}{"sv"};
			$svindex{$k}->[$i] = $score;
		}
	}
	$i++;
	%gh = ();
}

#print Dumper(\%sindex);

# score
foreach my $k (keys %sindex){
	my $str_tmp = join "",@{$sindex{$k}};
	if($str_tmp =~ /NA/){
		next;
	}else{
		my $arr = $sindex{$k};
		my @oarr;
		my @msdarr;
		for(my $i = 0; $i < scalar(@$arr); $i++){
			my @tarr;
			for(my $j = 0; $j < scalar(@$arr); $j++){
				if($i != $j){
					push @tarr,$arr->[$j];
				}
			}
			my @r;
			calOneLeaveOutZscore(\@tarr,\@r);
			push @msdarr,@r;
			my $ozscore = ($arr->[$i] - $r[0])/$r[1];
			push @oarr,$ozscore;
		}
	
		# print out mean and sd	
		#print "$k\t";
		#print join "\t",@oarr;
		##print "\t";
		##print join "\t",@$arr;
		#print "\n";
		
		# print out leave one out zscore
		my $flag_count = 0;
		foreach my $z (@oarr){
			if($z > 1.96 || $z < -1.96){
				$flag_count++;
			}
		}
		
		if($flag_count == 0){
		## output one leave out zscore
		#print "$k\t";
		#print join "\t",@oarr;
		#print "\n";\
		## output mean and sd for case comparision
		print "$k\t";
		print join "\t",@msdarr;
		print "\n";
		}
	}
}

sub calOneLeaveOutZscore {
	my $larr = shift;
	my $zarr = shift;
	my $lstd = std_dev(@$larr);
	my $lmean = average(@$larr);
	push @$zarr,$lmean;
	push @$zarr,$lstd;
	#my $sample_count = scalar(@$larr);
	#my $lmean = sum(@$larr)/$sample_count;
	
}

sub average {
        my (@values) = @_;

        my $count = scalar @values;
        my $total = 0; 
        $total += $_ for @values; 

        return $count ? $total / $count : 0;
}

sub std_dev {
        my ($average, @values) = @_;

        my $count = scalar @values;
        my $std_dev_sum = 0;
        $std_dev_sum += ($_ - $average) ** 2 for @values;

        return $count ? sqrt($std_dev_sum / $count) : 0;
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
