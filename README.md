# Urban Rural Pollen

This archive contains the bioinformatic commands, reference sequence databases and supplementary materials associated with the manuscript, 'Honey bee pollen foraging ecology across an urbanization gradient.' Refer to the Workflow.md file for an overview of the analysis and step by step commands. The empirical amplicon sequence data produced in this work are available at NCBI SRA (accession: PRJNA574766). The OhioDatabases.zip file contains reference sequence sets and associated taxonomic lineages for the _trnL_, _rbcL_ and ITS2 markers. These databases were used to assign taxonomy to the amplicon sequence data. Sequences present in these databases were downloaded from NCBI Nucleotide and filtered to only include sequences from species present in Ohio according to the USDA database. These geographically filtered sequence sets are part of a larger curation effort of the global corpus of these marker sequence data, described in detail in http://dx.doi.org/10.1101/672782 and available at https://github.com/RTRichar/MetabarcodeDBsV2. 

### Overview of contents
1. Sample-tag indices and PCR primers  
- Coming soon

2. Bioinformatic analysis  
-Workflow.md contains overview of tools and commands used to process sequencing data
-OhioDatabases.zip contain reference sequence databases used for taxonomic annotationo of sequence data
-Python code used to format and summarize taxonomic annotations: VsearchToMetaxa2.py and ThreeMarkerMedianConsensus.py

3. Directories containing Excel spreadsheets of raw data, including genus and family-leel three-marker metabarcoding median results as well as rarefied ITS2 results  
-Family-leel three-marker metabarcoding median results: FamilyMedians
-Genus-leel three-marker metabarcoding median results: GenusMedians
-Genus-level rarefied ITS2 results: RarefiedITS2
