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

        Select regions with differ fragment length diversity (measured by PFE) in 2 groups.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 s1.pfe,s2.pfe...[group1] s3.pfe,s5.pfe...[group2]

USAGE
print "$usage";
exit(1);
};

print STDERR `date`;
print STDERR "analysis start\n";

my $g1_str = shift;
my $g2_str = shift;

my @g1 = split /,/,$g1_str;
my @g2 = split /,/,$g2_str;

# build the index
my @g = (@g1, @g2);

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

#print Dumper(\%svindex);
#
# selection
# vscore
##foreach my $k (keys %svindex){
##	my $str_tmp = join "",@{$svindex{$k}};
##	if($str_tmp =~ /NA/){
##		next;
##	}else{
##		my $avg = abs(sum(@{$svindex{$k}})/scalar(@{$svindex{$k}}));
##		if($avg < 0.2){
##			print $k,"\t";
##			print join "\t",@{$svindex{$k}};
##			print "\n";
##		}
##	}
##}

# score
foreach my $k (keys %sindex){
	my $str_tmp = join "",@{$sindex{$k}};
	if($str_tmp =~ /NA/){
		next;
	}else{
		my $arr = $sindex{$k};
		my $n_avg = (abs($arr->[0]) + abs($arr->[1]) + abs($arr->[2]))/3;
		my $t_avg = (abs($arr->[3]) + abs($arr->[4]) + abs($arr->[5]))/3;
		if($arr->[3] > 0 && $arr->[4] > 0 && $arr->[5] > 0 && ($t_avg+0.01)/($n_avg+0.01) > 2){
		#if($arr->[3] > 0 && $arr->[4] > 0 && $arr->[5] > 0 && $n_avg >= 0.1 && $t_avg >= 0.1){
			my @k = split /#/,$k;
			print join "\t",@k;
			print "\t";
			print join "\t",@{$sindex{$k}};
			print "\n";
		}
	}
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
		$h->{$id}{"s"} = $score;
		$h->{$id}{"sv"} = abs($vscore);
		#print "$chr\t$s\t$e\t$n1\t$n2\t$strand\t$score\t$vscore\n";
	}
	close IN;
}

print STDERR "analysis done\n";
print STDERR `date`;
