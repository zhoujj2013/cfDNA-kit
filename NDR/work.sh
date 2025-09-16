perl -ne 'chomp; my @t = split /\t/; print "/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python ~/github/cfDNA-kit/NDR/ndr.py --bed_file ../ref/gencode.v45.annotation.tss.b1k.bed --bam_file $t[1] --o $t[0].ndr.tsv\n";' samples.lst > cal.ndr.sh

ls *.ndr.tsv | while read line; do echo "perl sum.pl $line > $(basename $line ".tsv").sum.tsv";done > sum.sh;

perl ~/bin/joinTable.pl PL230613001.ndr.sum.tsv PL230613002.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613003.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613004.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613005.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613006.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613007.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613008.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613009.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613010.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613011.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613012.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613013.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613014.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613015.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613016.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613017.ndr.sum.tsv -of 1.all,2.2 | perl ~/bin/joinTable.pl - PL230613018.ndr.sum.tsv -of 1.all,2.2 > ndr.sum.tsv

python ./ttest.py ndr.sum.tsv ndr.ttest_ind.tsv

less ndr.ttest_ind.tsv | awk '$3 < 0.01 && $5 < 0.01' | sort -k2n | awk '$2 < 0' | less -S
