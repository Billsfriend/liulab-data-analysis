for file in ~/Documents/ZCMiceBulk/sorted/*.sorted.bam; do Rscript featurecounts.R $file ~/Documents/ZCMiceBulk/ref/mouse.genome.gtf 10 $file; done

# merge counts
perl -lne 'if ($ARGV=~/(.*).count/){print "$1\t$_"}' *.count > merge.count

# extract counts value and tpm value
awk -F"\t" '{print $1"\t"$2"\t"$3}' merge.count > gene.count

awk -F"\t" '{print $1"\t"$2"\t"$5}' merge.count > gene.tpm


# https://www.jianshu.com/p/3352bfd404f3
