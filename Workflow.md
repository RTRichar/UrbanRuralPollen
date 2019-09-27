# 
#### for rbcL and ITS2
'''
pear -j 12 -y 48000M -f SampleName_R1.fastq -r SampleName_R2.fastq -o SampleName_PEARED.fastq -t 100 -q 20 -n 250
'''

#### for trnL only
'''
metaxa2 -1 SampleName_R1.fastq -2 SampleName_R2.fastq -o ${a}_trnL_MTXA -g OHIO_TRNL_G --taxlevel 6 --cpu 28 --distance 10
'''

