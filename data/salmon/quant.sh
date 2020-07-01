#!/bin/bash
cd rawdata/
for R1 in *R1_001.fastq.gz
do
R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
samp=$(basename $R1 _R1_001.fastq.gz)
echo $R1 $R2 $samp
echo "Processing sample ${samp}"
salmon quant -i salmon_index -l A \
         -1 ${R1} \
         -2 ${R2} \
         -p 8 --validateMappings -o ../quants/${samp}_quant
done 
