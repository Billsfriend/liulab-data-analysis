#!/bin/bash
for fn in ~/Documents/DataAnalysis/ZCMiceBulk/data/B{1..3};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i salmon_mouse_index -l A \
         -1 ${fn}/${samp}_1.clean.fq.gz \
         -2 ${fn}/${samp}_2.clean.fq.gz \
         -p 8 --validateMappings --gcBias -o quants/${samp}_quant
done
