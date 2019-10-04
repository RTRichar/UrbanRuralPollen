#

## Concatenate forward and reverse reads into single contiguous sequence using Pear or Metaxa2
```
# for rbcL and ITS2
pear -j 12 -y 48000M -f SampleName_R1.fastq -r SampleName_R2.fastq -o SampleName_PEARED.fastq -t 100 -q 20 -n 250

# for trnL only (since trnL amplicons are often longer than 600 bp, merging with Pear is not feasible)
metaxa2 -1 SampleName_R1.fastq -2 SampleName_R2.fastq -o SampleName_trnL_MTXA -g OHIO_TRNL_G --cpu 28 --distance 10 --truncate F --save_raw T
# after running, Metaxa2 produces a SampleName_raw directory which contains a fasta of the concatenated reads. Use this for downstream VSEARCH
```

## Query concatenated reads against curated databases 
```
# fastq must be converted to fasta before VSEARCH alignment 
vsearch --fastq_filter SampleName_PEARED.fastq -fastaout SampleName_PEARED.fasta
# rbcL
vsearch --usearch_global SampleName_PEARED.fasta --db rbcLOhioGenera.fa --id 0.96 --maxaccepts 100 --maxrejects 50 --maxhits 1 --gapopen 0TE --gapext 0TE --userout SampleName_rbcL.txt --alnout SampleName_rbcL.aln --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov --threads 28 --target_cov 0.98
# ITS2
vsearch --usearch_global SampleName_PEARED.fasta --db ITS2OhioGenera.fa --id 0.925 --maxaccepts 100 --maxrejects 50 --maxhits 1 --gapopen 0TE --gapext 0TE --userout SampleName_ITS2.txt --alnout SampleName_ITS2.aln --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov --threads 28 --target_cov 0.98
# trnL
vsearch --usearch_global SampleName_Mtxa2_Raw.fasta --db trnLOhioGenera.fa --iddef 0 --maxaccepts 200 --maxrejects 200 --maxhits 1 --gapopen 0I/0TE --gapext 0TE --userout SampleName_trnL.txt --alnout SampleName_trnL.aln --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov --id 0.8 --threads 28 --target_cov 0.8
```

## Use python script to parse alignment files, assign sequences to taxonomic lineages and produce summary tables
```
VsearchToMetaxa2.py -v SampleName_trnL.txt -t trnLOhioGenera.tax -o SampleName.trnL.taxonomy.txt
VsearchToMetaxa2.py -v SampleName_rbcL.txt -t rbcLOhioGenera.tax -o SampleName.rbcL.taxonomy.txt
VsearchToMetaxa2.py -v SampleName_ITS2.txt -t ITS2OhioGenera.tax -o SampleName.ITS2.taxonomy.txt
```

## Use python script to get consensus-filtered median proportional abundances
#### More detailed description at https://github.com/RTRichar/QuantitativePollenMetabarcoding
```
# for help screen
python ThreeMarkerMedianConsensus.py --help

# requires file specifying base filenames, one for each sample
head -n 3 FSR14_BaseNames.txt
4F10_S57
4F11_S68
4F12_S78

# requires file specifying marker names 
-bash-4.2$ more ./../Markers.txt
ITS2
trnL
rbcL

# input filename convention must be SampleBaseName.MarkerName.taxonomy.txt (e.g. 4F10_S57.ITS2.taxonomy.txt)

# example of trnL analysis of FSR data at genus rank 
python ThreeMarkerMedianConsensus.py -b FSR14_BaseNames.txt -m Markers.txt -r 5 -pf 0.0005 -mt 0.001 -o FSR
```
