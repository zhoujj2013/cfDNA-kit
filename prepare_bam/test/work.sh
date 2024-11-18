gunzip -d /mnt/dfc_data1/home/lixuefei/Project/20200914TPcf/cfDNA12/RawData/211019_SEQ014_FP100002978BR_L01_PL210929002001/FP100002978BR_L01_38_2.fq.gz -c > FP100002978BR_L01_38_2.fq &
gunzip -d /mnt/dfc_data1/home/lixuefei/Project/20200914TPcf/cfDNA12/RawData/211019_SEQ014_FP100002978BR_L01_PL210929002001/FP100002978BR_L01_38_1.fq.gz -c > FP100002978BR_L01_38_1.fq &

perl ../bin/cfdna_mNGS.pl config.txt
