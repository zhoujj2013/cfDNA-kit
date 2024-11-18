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

        This script is designed for mNGS analysis (species identification).
        Author: zhoujj2013\@gmail.com 
        Usage: $0 <config.txt>

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;

my $sample_f = abs_path($conf{SAMPLE});
my $out = abs_path($conf{OUTDIR});
my $thread = $conf{THREAD};

#### write you things ###
#decide to add function in another file
#my %compare;
open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^#/);
	my @t = split /\t/;
	my $id = shift @t;

	# determine PE
	my $type = "PE";
	my ($r1,$r2);
	if(scalar(@t) == 1){
		$type ="SE";
		$r1 = abs_path($t[0]);
		$r2 = $r1;
	}elsif(scalar(@t) == 2){
		$type = "PE";
		$r1 = abs_path($t[0]);
		$r2 = abs_path($t[1]);
	}
		
	## QC start
	make_path "$out/$id/00datafilter";
	my $encode_type = `perl $conf{BIN}/phredDetector.pl $r1`;
	my $filter_cutoff = 5;
	my $phred = "";
    	if($encode_type == 33){
        	#$filter_cutoff = 38;
        	$phred = " ";
    	}elsif($encode_type == 64){
        	#$filter_cutoff = 69;
        	$phred = " --phred64 ";
    	}
	if($type eq "PE"){
	$mk .= "00fastqc.origin.finished: $r1 $r2\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $r1 $r2 > /dev/null 2>/dev/null && touch 00fastqc.origin.finished\n";
	$all .= "00fastqc.origin.finished ";

	$mk .= "00trim.finished: 00fastqc.origin.finished\n";	
	$mk .= "\t$conf{AdapterRemoval} --qualitybase $encode_type  --file1 $r1 --file2 $r2 --settings  $out/$id/00datafilter/$id.adapterRemoval.setting --output1 $out/$id/00datafilter/$id.Trim.R1.fq --output2 $out/$id/00datafilter/$id.Trim.R2.fq --singleton $out/$id/00datafilter/$id.Single.fq --discarded $out/$id/00datafilter/$id.Discarded.fq --minlength $conf{MINLEN} --trimns --trimqualities --minquality $filter_cutoff --trimwindows 15 --threads $conf{THREAD} --adapter-list $conf{ADAPTER} > $out/$id/00datafilter/$id.trim.log 2>$out/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	$all .= "00trim.finished ";

	$mk .= "00fastqc.finished: 00trim.finished\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $out/$id/00datafilter/$id.Trim.R1.fq $out/$id/00datafilter/$id.Trim.R2.fq > /dev/null 2>/dev/null && touch 00fastqc.finished\n";
	$all .= "00fastqc.finished ";

	}elsif($type eq "SE"){

	$mk .= "00fastqc.origin.finished: $r1\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $r1 > /dev/null 2>/dev/null && touch 00fastqc.origin.finished\n";
	$all .= "00fastqc.origin.finished ";

	$mk .= "00trim.finished: 00fastqc.origin.finished\n";	
	$mk .= "\t$conf{AdapterRemoval} --qualitybase $encode_type  --file1 $r1 --settings  $out/$id/00datafilter/$id.adapterRemoval.setting --output1 $out/$id/00datafilter/$id.Trim.R1.fq --singleton $out/$id/00datafilter/$id.Single.fq --discarded $out/$id/00datafilter/$id.Discarded.fq --minlength $conf{MINLEN} --trimns --trimqualities --minquality $filter_cutoff --trimwindows 15 --threads $conf{THREAD} --adapter-list $conf{ADAPTER} > $out/$id/00datafilter/$id.trim.log 2>$out/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	$all .= "00trim.finished ";

	$mk .= "00fastqc.finished: 00trim.finished\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $out/$id/00datafilter/$id.Trim.R1.fq > /dev/null 2>/dev/null && touch 00fastqc.finished\n";
	$all .= "00fastqc.finished ";
	}
	## QC end ##
		
	## alignment
	my $single = "$out/$id/00datafilter/$id.Single.fq";
	my $discarded = "$out/$id/00datafilter/$id.Discarded.fq";
	my $clear_r1 = "$out/$id/00datafilter/$id.Trim.R1.fq";
	my $clear_r2 = "$out/$id/00datafilter/$id.Trim.R2.fq";
	
	my $reads = "";
	if($type eq "PE"){
		$reads = "$clear_r1 $clear_r2";
	}elsif($type eq "SE"){
		$reads = "$clear_r1";
	}
	
	make_path "$out/$id/01alignment";
	$mk .= "01align.finished: 00fastqc.finished\n";
	# bwa alignment
	$mk .= "\t$conf{BWA} mem -t $thread -k 36 -M $conf{INDEX} $reads 2>/dev/null | $conf{SAMTOOLS} view -Sb -@ $thread -T $conf{INDEX} -o $out/$id/01alignment/$id.BackgroudAlign.bam > /dev/null 2>/dev/null ";
	$mk .= "&& touch 01align.finished\n";
	$all .= "01align.finished ";

	# retrive unmapped reads (relative to ref.)
	make_path "$out/$id/01alignment/TMP" unless(-d "$out/$id/01alignment/TMP");
	$mk .= "01getUnmapped.finished: 01align.finished\n";
	$mk .= "\t$conf{GATK} SortSam --INPUT $out/$id/01alignment/$id.BackgroudAlign.bam --OUTPUT $out/$id/01alignment/$id.BackgroudAlign.sorted.bam --SORT_ORDER coordinate --TMP_DIR $out/$id/01alignment/TMP > /dev/null 2>/dev/null ";
	$mk .= "&& $conf{SAMTOOLS} index $out/$id/01alignment/$id.BackgroudAlign.sorted.bam > /dev/null 2>/dev/null ";
	$mk .= "&& $conf{SAMTOOLS} flagstat $out/$id/01alignment/$id.BackgroudAlign.sorted.bam > $out/$id/01alignment/$id.BackgroudAlign.sorted.bam.flagstat ";
	$mk .= "&& $conf{SAMTOOLS} view -b -f 4 $out/$id/01alignment/$id.BackgroudAlign.sorted.bam > $out/$id/01alignment/$id.BackgroudAlign.sorted.unmapped.bam ";

	if($type eq "PE"){
	$mk .= "&& $conf{BEDTOOLS} bamtofastq -i $out/$id/01alignment/$id.BackgroudAlign.sorted.unmapped.bam -fq $out/$id/01alignment/$id.unmappedR1.fastq -fq2 $out/$id/01alignment/$id.unmappedR2.fastq > /dev/null 2>/dev/null ";
	}elsif($type eq "SE"){
	$mk .= "&& $conf{BEDTOOLS} bamtofastq -i $out/$id/01alignment/$id.BackgroudAlign.sorted.unmapped.bam -fq $out/$id/01alignment/$id.unmappedR1.fastq  > /dev/null 2>/dev/null ";
	}

	$mk .= "&& touch 01getUnmapped.finished\n";
	$all .= "01getUnmapped.finished ";

	# get insert and duplication rate
	$mk .= "01getInsertAndDup.finished: 01getUnmapped.finished\n";
	$mk .= "\tjava -jar $conf{PICARD} CollectInsertSizeMetrics INPUT=$out/$id/01alignment/$id.BackgroudAlign.sorted.bam OUTPUT=$out/$id/01alignment/$id.insert.txt HISTOGRAM_FILE=$out/$id/01alignment/$id.insert.hist.pdf >$out/$id/01alignment/insert.log 2>$out/$id/01alignment/insert.err ";
	$mk .= "&& java -jar $conf{PICARD} MarkDuplicates I=$out/$id/01alignment/$id.BackgroudAlign.sorted.bam O=$out/$id/01alignment/$id.sorted.rmdup.bam M=$out/$id/01alignment/$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ >$out/$id/01alignment/rmdup.log 2>$out/$id/01alignment/rmdup.err ";
	$mk .= "&& touch 01getInsertAndDup.finished\n";
	$all .= "01getInsertAndDup.finished ";
	
	make_path "$out/$id/02kraken";
	my $is_paired = "--paired";
	my $depleted_host = "$out/$id/01alignment/$id.unmappedR1.fastq $out/$id/01alignment/$id.unmappedR2.fastq";
	if($type eq "SE"){
		$is_paired = "";
		$depleted_host = "$out/$id/01alignment/$id.unmappedR1.fastq";
	}
	# database download
	# https://benlangmead.github.io/aws-indexes/k2
	$mk .= "02kraken.finished: 01getUnmapped.finished\n";
	$mk .= "\texport KRAKEN2_DB_PATH=\"$conf{KRAKEN2_DB_PATH}:\" && $conf{kraken2} --db $conf{kdb} --threads $thread $is_paired --classified-out $out/$id/02kraken/$id.cseqs#.fq --unclassified-out $out/$id/02kraken/$id.uncseqs#.fq --use-names $depleted_host > $out/$id/02kraken/$id.spe.count 2> $out/$id/02kraken/$id.kraken.log ";
	$mk .= "&& cut -f 3 $out/$id/02kraken/$id.spe.count | sort | uniq -c | sort -k1nr > $out/$id/02kraken/$id.spe.count.stat && perl $conf{BIN}/convert_whitespace2tab.pl $out/$id/02kraken/$id.spe.count.stat > $out/$id/02kraken/$id.spe.count.stat.rf ";
	$mk .= "&& touch 02kraken.finished\n";
	$all .= "02kraken.finished ";
	## kraken --- end

	make_path abs_path($conf{OUTDIR});
	open OUT, ">$out/$id/makefile";
	print OUT $all, "\n";
	print OUT $mk, "\n";
	close OUT;
	$all = "all: ";
	$mk = "";
}
close IN;

#########################
sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
