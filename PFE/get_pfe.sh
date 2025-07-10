perl ~/bin/fishInWinter.pl -fc 4 ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/Housekeeping/HK_gene_names.txt ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > HK_gene_names.bed

perl ~/gitee/ngskit/ngsplot/getTss8bed.pl HK_gene_names.bed > HK_gene_names.tss.bed

############
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK_gene_names.tss.bed 0 1000 > HK_gene_names.0.1000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK_gene_names.tss.bed -1000 0 > HK_gene_names.m1000.0.pfe

python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK_gene_names.tss.bed -5000 -4000 > HK_gene_names.m5000.m4000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK_gene_names.tss.bed 4000 5000 > HK_gene_names.4000.5000.pfe


#############
#cat ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/Housekeeping/HK_gene_names.txt ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/FANTOM5/Fantom5_all_lower0.1.txt | perl ~/bin/fishInWinter.pl -fc 4 - ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > HK.lowerE.genes.bed

cat ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/FANTOM5/Fantom5_all_lower0.1.txt | perl ~/bin/fishInWinter.pl -fc 4 - ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > HK.lowerE.genes.bed

perl ~/gitee/ngskit/ngsplot/getTss8bed.pl HK.lowerE.genes.bed > HK.lowerE.genes.tss.bed


python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK.lowerE.genes.tss.bed 0 1000 > HK.lowerE.0.1000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK.lowerE.genes.tss.bed -1000 0 > HK.lowerE.m1000.0.pfe

python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK.lowerE.genes.tss.bed -5000 -4000 > HK.lowerE.m5000.m4000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam HK.lowerE.genes.tss.bed 4000 5000 > HK.lowerE.4000.5000.pfe

##############
perl ~/gitee/ngskit/ngsplot/getTss8bed.pl ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > all.genes.tss.bed
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam all.genes.tss.bed 0 1000 > all.genes.0.1000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam all.genes.tss.bed -1000 0 > all.genes.m1000.0.pfe

python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam all.genes.tss.bed -5000 -4000 > all.genes.m5000.m4000.pfe
python3 ../../analysis/cal_pfe.py bam221/HC_zhou2.BackgroudAlign.sorted.bam all.genes.tss.bed 4000 5000 > all.genes.4000.5000.pfe

perl ~/bin/joinTable.pl -k1 4 -k2 4 HK.lowerE.m5000.m4000.pfe HK.lowerE.m1000.0.pfe -of 1.all,2.7 | perl ~/bin/joinTable.pl -k1 4 -k2 4 - HK.lowerE.0.1000.pfe  -of 1.all,2.7 | perl ~/bin/joinTable.pl -k1 4 -k2 4 - HK.lowerE.4000.5000.pfe -of 1.all,2.7 | grep -v "nan" | grep -v "none" > HK.lowerE.pfe

perl ~/bin/joinTable.pl -k1 4 -k2 4 HK_gene_names.m5000.m4000.pfe HK_gene_names.m1000.0.pfe -of 1.all,2.7 | perl ~/bin/joinTable.pl -k1 4 -k2 4 - HK_gene_names.0.1000.pfe  -of 1.all,2.7 | perl ~/bin/joinTable.pl -k1 4 -k2 4 - HK_gene_names.4000.5000.pfe -of 1.all,2.7 | grep -v "nan" | grep -v "none" > HK_gene_names.pfe

Rscript  ~/github/jjUtil/compareDisRNAseq/MultiBoxplotNoOutlineSize.r HK_gene_names 0.5 6 6 HK_gene_names.pfe:7:darkgray HK_gene_names.pfe:8:darkgray HK_gene_names.pfe:9:darkgray HK_gene_names.pfe:10:darkgray

Rscript  ~/github/jjUtil/compareDisRNAseq/MultiBoxplotNoOutlineSize.r HKlowerE 0.5 6 6 HK.lowerE.pfe:7:darkgray HK.lowerE.pfe:8:darkgray HK.lowerE.pfe:9:darkgray HK.lowerE.pfe:10:darkgray

Rscript  ~/github/jjUtil/compareDisRNAseq/MultiBoxplotNoOutlineSize.r HKlowerE.vs.HKhighE 0.5 6 6 HK.lowerE.pfe:7:orangered HK_gene_names.pfe:7:blue HK.lowerE.pfe:8:orangered HK_gene_names.pfe:8:blue HK.lowerE.pfe:9:orangered HK_gene_names.pfe:9:blue HK.lowerE.pfe:10:orangered HK_gene_names.pfe:10:blue

