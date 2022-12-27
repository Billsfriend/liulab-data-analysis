genome= ~/Documents/ZCMiceBulk/ref/mouse.genome.fna

genome_prefix= ~/Documents/ZCMiceBulk/ref/mouse

sudo hisat2-build $genome $genome_prefix 1>hisat2-build.log 2>hisat2-build.err

