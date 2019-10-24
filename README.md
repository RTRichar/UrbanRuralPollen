# Urban Rural Pollen

### About this archive
This archive contains the bioinformatic commands, reference sequence databases and supplementary materials associated with the manuscript, 'Honey bee pollen foraging ecology across an urbanization gradient.' Refer to the Workflow.md file for an overview of the analysis and step by step commands. The empirical amplicon sequence data produced in this work are available at NCBI SRA (accession: PRJNA574766). The OhioDatabases.zip file contains reference sequence sets and associated taxonomic lineages for the _trnL_, _rbcL_ and ITS2 markers. These databases were used to assign taxonomy to the amplicon sequence data. Sequences present in these databases were downloaded from NCBI Nucleotide and filtered to only include sequences from species present in Ohio according to the USDA database. These geographically filtered sequence sets are part of a larger curation effort of the global corpus of these marker sequence data, described in detail in http://dx.doi.org/10.1101/672782 and available at https://github.com/RTRichar/MetabarcodeDBsV2. 

### More detailed overview of contents
1. Sample-tag indices and PCR primers  
    - __SampleTag_SiteDate_Index.xlsx__ contains information on the relationship between NCBI file label sample tags and sample site date
    - __PCR_methods.xlsx__ contains information on the primers and PCR conditions used for library preparation

2. Bioinformatic analysis  
    - __Workflow.md__ contains overview of tools and commands used to process sequencing data
    - __OhioDatabases.zip__ contains reference sequence databases used for taxonomic annotation of sequence data
    - __VsearchToMetaxa2.py__ and __ThreeMarkerMedianConsensus.py__ were used to format and summarize taxonomic annotations

3. Directories containing Excel spreadsheets of raw data  
     - Family-leel three-marker metabarcoding median results can be found in the __FamilyMedians__ directory
     - Genus-leel three-marker metabarcoding median results can be found in the __GenusMedians__ directory
     - Genus-level rarefied ITS2 results can be found in the __RarefiedITS2__ directory
