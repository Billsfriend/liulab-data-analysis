for file in *.sam; do 

    samtools view -@ 10 -b $file > ~/Documents/ZCMiceBulk/sorted/${file%.sam}.bam

    samtools sort -@ 10 ~/Documents/ZCMiceBulk/sorted/${file%.sam}.bam > ~/Documents/ZCMiceBulk/sorted/${file%.sam}.sorted.bam

    samtools index -@ 10 ~/Documents/ZCMiceBulk/sorted/${file%.sam}.sorted.bam

done

# https://www.jianshu.com/p/3352bfd404f3
