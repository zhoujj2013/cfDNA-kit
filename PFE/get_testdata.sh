perl ~/bin/fishInWinter.pl -fc 4 ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/Housekeeping/HK_gene_names.txt ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > HK_gene_names.bed

perl ~/gitee/ngskit/ngsplot/getTss8bed.pl HK_gene_names.bed > HK_gene_names.tss.bed

cat ../../data/blood_expressed/Nucleosome_ctDNA-master/ref/FANTOM5/Fantom5_all_lower0.1.txt | perl ~/bin/fishInWinter.pl -fc 4 - ~/project/01lncrna_structure/10NEI/hg38_v24/gencode.v24.annotation.bed > HK.lowerE.genes.bed

perl ~/gitee/ngskit/ngsplot/getTss8bed.pl HK.lowerE.genes.bed > HK.lowerE.genes.tss.bed
