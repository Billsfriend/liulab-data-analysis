genome_prefix=~/Documents/ZCMiceBulk/ref/mouse

# 写一个for循环进行批量运算

#for paired -end, PE

for sample in `perl -lpe 's/_1.clean.fq.gz//'`; do hisat2 --new-summary -p 2 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ${sample}_1.clean.fq.gz -2 ${sample}_2.clean.fq.gz -S ${sample}.sam 2> ${sample}.err; done

# do it one by one
hisat2 --new-summary -p 6 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ~/Documents/ZCMiceBulk/data/ZC-A2_FKDL220090431-1a-AK438_1.clean.fq.gz -2 ~/Documents/ZCMiceBulk/data/ZC-A2_FKDL220090431-1a-AK438_2.clean.fq.gz -S ~/Documents/ZCMiceBulk/hisat2/ZC-A2.sam 2> ~/Documents/ZCMiceBulk/hisat2/ZC-A2.err

# -p 8 means run in 8 threads parallely. CPU 80% occupied.
hisat2 --new-summary -p 8 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ~/Documents/ZCMiceBulk/data/ZC-A3_FKDL220090431-1a-AK325_1.clean.fq.gz -2 ~/Documents/ZCMiceBulk/data/ZC-A3_FKDL220090431-1a-AK325_2.clean.fq.gz -S ~/Documents/ZCMiceBulk/hisat2/ZC-A3.sam 2> ~/Documents/ZCMiceBulk/hisat2/ZC-A3.err

hisat2 --new-summary -p 10 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ~/Documents/ZCMiceBulk/data/ZC-B1_FKDL220090431-1a-AK434_1.clean.fq.gz -2 ~/Documents/ZCMiceBulk/data/ZC-B1_FKDL220090431-1a-AK434_2.clean.fq.gz -S ~/Documents/ZCMiceBulk/hisat2/ZC-B1.sam 2> ~/Documents/ZCMiceBulk/hisat2/ZC-B1.err

hisat2 --new-summary -p 10 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ~/Documents/ZCMiceBulk/data/ZC-B2_FKDL220090431-1a-AK435_1.clean.fq.gz -2 ~/Documents/ZCMiceBulk/data/ZC-B2_FKDL220090431-1a-AK435_2.clean.fq.gz -S ~/Documents/ZCMiceBulk/hisat2/ZC-B2.sam 2> ~/Documents/ZCMiceBulk/hisat2/ZC-B2.err

hisat2 --new-summary -p 10 -x ~/Documents/ZCMiceBulk/ref/mouse -1 ~/Documents/ZCMiceBulk/data/ZC-B3_FKDL220090431-1a-AK436_1.clean.fq.gz -2 ~/Documents/ZCMiceBulk/data/ZC-B3_FKDL220090431-1a-AK436_2.clean.fq.gz -S ~/Documents/ZCMiceBulk/hisat2/ZC-B3.sam 2> ~/Documents/ZCMiceBulk/hisat2/ZC-B3.err