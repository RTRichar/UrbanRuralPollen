#!/usr/bin/env python
'''
--$ python Analyze_Multi-Locus.py -b base_names.txt -m marker_names.txt -pf pre-filter_threshold -mt masking_threshold -o output_name.csv

This script takes the following inputs: (1) an input file containing the base names for each sample; (2) a file of the markers used, one marker per line; (3) a protortional threshold, between 0 and 1, indicating the minimum proportional threshold for a family to be considered present in the sample; (4) a protortional threshold, between 0 and 1, indicating the median proportional abundance a taxonomic group needs to be represented by in order to be included in the final output (e.g. 0.05 means that any taxon below 5% abundance will be excluded); (5) output file name

This script has to be run in the same directory as MTXA2 output files and all files need to have the same naming convention (sampleID.MarkerName.taxonomy.txt, with sampleIDs being identical to those provided in sys.argv[1] and MarkerName being identical to those provided in sys.argv[2])
'''

from __future__ import division # for decimals, not sure if needed
import re # for regular expressions
import sys # standard input/output
import argparse
import numpy
def median(lst):
    return numpy.median(numpy.array(lst)) # for median calculation

parser = argparse.ArgumentParser()
# required
parser.add_argument('-b', '--BaseFileNames', required = True, help = "File containing base file names of all sample files")
parser.add_argument('-m', '--MarkerNames', required = True, help = "File containing name o each marker") 
parser.add_argument('-o', '--OutPutPrefix', required = True, help = "Prefix name of output csv files")
parser.add_argument('-r', '--Rank', required = True, help = "Rank to be analyzed (INT from 0 to 6 [kingdom to species])")
# optional
parser.add_argument('-pf', '--PreFilterThreshold', required = True, help = "Proportion abundance of sequences required for a taxon to be included in the metabarcoding median calculation")
parser.add_argument('-mt', '--MaskThreshold', required = True, help = "Proportion abundance of sequences required for a taxon to be written to output csv file")
args = parser.parse_args()

# make list to define what rank is to be analyzed
Ranks = ['King','Phyl','Class','Order','Fam','Gen','Sp']

ifile = open(args.BaseFileNames, 'r') # base names
mfile = open(args.MarkerNames, 'r') # markers
pfilt_thresh = float(args.PreFilterThreshold)
msk_thresh = float(args.MaskThreshold) # Taxanomic groups are masked into 'other' category if their abundance is below this threshold
MedOutName = str(args.OutPutPrefix + '_Med_' + Ranks[int(args.Rank)] + 'Pf' + args.PreFilterThreshold + 'Mt' + args.MaskThreshold + '.csv')
#print pfilt_thresh, msk_thresh

markers = mfile.read().splitlines() # make list of markers used in study
samples = ifile.read().splitlines() # make list of base names of samples
#print markers
#print samples

### Section 1.) make a dictionary with each key being the sample name and each value being a list of marker files for the sample key
Sample_Sets_Dct = {}
for a in samples: # for each base name and each marker, create filename and family dictionary name and add them to above lists
	a = str(a.split('.')[:1])[2:-2] 
	Set_Lst = []
	for b in markers: # iterate over each marker
		SampleMarker = str(a + '.' + b + '.taxonomy.txt') # variable containing file name
		Set_Lst.append(SampleMarker)
	Sample_Sets_Dct[a] = Set_Lst

### Section 2.) make four dictionaries, one containing lists of count dictionaries per sample key, one containing lists of proportion dictionaries per sample key, one containing lists of median proportion dictionaries per sample key for families detected in each sample, and one containing lists of median proportion dictionaries per sample key for families detected accross all samples
File_Count_Dct = {}
File_Prop_Dct = {}
File_Med_Dct = {}
File_Med_Msk_Dct = {}
File_Med_AllFams_Dct = {}
for c in samples: # for each base name and each marker, create filename and family dictionary name and add them to above lists
	c = str(c.split('.')[:1])[2:-2]
	Count_Set_Lst = []
	Prop_Set_Lst = []
	for d in markers: # iterate over each marker
		SampleFileCount = str(c + '_' + d + '_count') # variable containing file name
		SampleFileProp = str(c + '_' + d + '_prop')
		Count_Set_Lst.append(SampleFileCount)
		Prop_Set_Lst.append(SampleFileProp)
	File_Count_Dct[c] = Count_Set_Lst
	File_Prop_Dct[c] = Prop_Set_Lst
for c in samples:
        c = str(c.split('.')[:1])[2:-2]
        SampleFileMedMsk = str(c + '_medMsk')
        File_Med_Msk_Dct[c] = SampleFileMedMsk
for c in samples:
        c = str(c.split('.')[:1])[2:-2]
        SampleFileMed = str(c + '_med')
        File_Med_Dct[c] = SampleFileMed
for c in samples:
        c = str(c.split('.')[:1])[2:-2]
        SampleFileMedAllFam = str(c + '_med_AllFams')
        File_Med_AllFams_Dct[c] = SampleFileMedAllFam

### Section 3.) populate count dictionaries
for d in samples:
	for e in range(3):
		File_Count_Dct[d][e] = {}
		with open(Sample_Sets_Dct[d][e], 'r') as MTXA_file:
			for line in MTXA_file:
				line = line.split('\t')[1]
				if re.search( r'k__', line) is not None: # read only the lines containing an assignment
					if len(line.split(';')) >= int(int(args.Rank) + int(1)):
						fam = str(line.split(';')[(int(args.Rank) - 1)] + '_' + line.split(';')[int(args.Rank)]) # add sub-rank	
						if fam not in File_Count_Dct[d][e]:
							File_Count_Dct[d][e][fam] = 1
						else:
							File_Count_Dct[d][e][fam] += 1 # add counts of each family to family counts dictionary for each input file
		MTXA_file.close()

### Section 4.) convert count dictionaries to prop dictionaties
for f in samples:
	for g in range(3):
		File_Prop_Dct[f][g] = {}
		total_count = sum(File_Count_Dct[f][g].values())
		for key, value in File_Count_Dct[f][g].iteritems():
			if float(value)/total_count > pfilt_thresh:
				File_Prop_Dct[f][g][key] = float(value)/total_count

### Section 4.1.)write data for each marker to file
for g in range(3):
        output_dictionary = {}
        families = []
	IndOut = str(args.OutPutPrefix + '_' + markers[g] + '_' + Ranks[int(args.Rank)] + 'Mt' + args.MaskThreshold + '.csv')
	markeroutfile = open(IndOut, 'w')
        header = []
        for f in samples:
                for key in File_Prop_Dct[f][g]:
                        if key not in families:
                                families.append(key)
                header.append(f)
        for fam in families:
                Lst = []
                for f in samples:
                        if fam in File_Prop_Dct[f][g]:
				if File_Prop_Dct[f][g][fam] >= float(msk_thresh):
	                                Lst.append(File_Prop_Dct[f][g][fam])
				else:
					Lst.append(0)
                        else:
                                Lst.append(0)
                output_dictionary[fam] = Lst
        output_dictionary['aaa_header'] = header
        for key, value in sorted(output_dictionary.items()):
                markeroutfile.write( str(key) + ',' + str(value).strip('[ | ]') + '\n' )
        markeroutfile.close()

### Section 5.) for each sample, find consensus families (families detected with at least two of the four barcodes)
for h in samples: # for each sample 
	File_Med_Dct[h] = {}
	ITS2_fams = list(File_Prop_Dct[h][0].keys())
	rbcL_fams = list(File_Prop_Dct[h][1].keys())
#	trnH_fams = list(File_Prop_Dct[h][2].keys())
	trnL_fams = list(File_Prop_Dct[h][2].keys()) # create list of families for each marker
	TMP_Consensus = [] # create temporary list to hold consensus families for this sample
	for i in ITS2_fams: # look at all pairwise permutations of the markers to find consensus families
		if i in rbcL_fams:
			TMP_Consensus.append(i)
#		elif i in trnH_fams:
#			if i not in TMP_Consensus:
#				TMP_Consensus.append(i)
		elif i in trnL_fams:
			if i not in TMP_Consensus:
				TMP_Consensus.append(i) # take ITS2 families and see if any are also detected in rbcL trnL or trnH. If yes, place family into 							       consensus list
	for j in rbcL_fams:
#		if j in trnH_fams:
#			if j not in TMP_Consensus:
#				TMP_Consensus.append(j)
		if j in trnL_fams:
			if j not in TMP_Consensus:
				TMP_Consensus.append(j) # take rbcL familes and see if any are also detected in trnH or trnL 
#	for k in trnH_fams:
#		if k in trnL_fams:
#			if k not in TMP_Consensus:
#				TMP_Consensus.append(k) # take trnH families and see if any are also detected in trnL
# describe
	for m in TMP_Consensus: # 1 2 4
		if m not in File_Med_Dct[h]:
			if m in File_Prop_Dct[h][0] and m in File_Prop_Dct[h][1] and m in File_Prop_Dct[h][2]:
				med = median([File_Prop_Dct[h][0][m], File_Prop_Dct[h][1][m], File_Prop_Dct[h][2][m]])
				File_Med_Dct[h][m] = med
	for m in TMP_Consensus: # 1 2
		if m not in File_Med_Dct[h]:
			if m in File_Prop_Dct[h][0] and m in File_Prop_Dct[h][1]:
				med = median([File_Prop_Dct[h][0][m], File_Prop_Dct[h][1][m]])
				File_Med_Dct[h][m] = med
	for m in TMP_Consensus: # 1 3
		if m not in File_Med_Dct[h]:
			if m in File_Prop_Dct[h][0] and m in File_Prop_Dct[h][2]:
				med = median([File_Prop_Dct[h][0][m], File_Prop_Dct[h][2][m]])
				File_Med_Dct[h][m] = med
	for m in TMP_Consensus: # 1 4
		if m not in File_Med_Dct[h]:
			if m in File_Prop_Dct[h][1] and m in File_Prop_Dct[h][2]:
				med = median([File_Prop_Dct[h][1][m], File_Prop_Dct[h][2][m]])
				File_Med_Dct[h][m] = med
# FOOTNOTE 1) the order of each four loop is important here because we have to get the medians of families found accross the samples according to how many samples they are found in (e.g. if they are found in four samples, we aren't interested in the median of only three samples)
# FOOTNOTE 2) for each sample and each consensus family,if family is detected with markers 1 through 4, take the median accross markers 1 through 4
# FOOTNOTE 3) for each sample and each consensus family,if family is detected with markers 1 through 3, take the median accross markers 1 through 3
# and so on and so forth for all pairwise permutations

### Section 6.) normalize median proportion dictionaries to sum of 1.0
for z in samples:
	prop_sum = sum(File_Med_Dct[z].values())
	if float(prop_sum) > 0:
		prop_dif = 1/float(prop_sum)
	for key, value in File_Med_Dct[z].iteritems():
			File_Med_Dct[z][key] = prop_dif*value

### Section 6.1.) mask all families below proportion threshold specified in sys.argv[3]
for z in samples:
        File_Med_Msk_Dct[z] = {}
        File_Med_Msk_Dct[z]['t__other'] = 0
        for key in File_Med_Dct[z]:
                if File_Med_Dct[z][key] > msk_thresh:
                        File_Med_Msk_Dct[z][key] = File_Med_Dct[z][key]
                else:
                        File_Med_Msk_Dct[z]['t__other'] += float(File_Med_Dct[z][key])

### Section 7.) find all consensus families accross all files and use them to populate all family medians dictionary
All_families = []
for z in samples:
	for key in File_Med_Msk_Dct[z]:
		if key not in All_families:
			All_families.append(key)
for x in samples:
	File_Med_AllFams_Dct[x] = {}
	for y in All_families:
		if y not in File_Med_Msk_Dct[x]:
			File_Med_AllFams_Dct[x][y] = 0
		else:
			File_Med_AllFams_Dct[x][y] = File_Med_Msk_Dct[x][y]
	
### Section 8.) use all family medians dictionary to populate final dictionary where every key is a family and evey value is a list of the median proportional abundances of each sample
Final_Dct = {}
for x in range(len(All_families)):
	name = All_families[x]
	All_families[x] = []
	for y in samples:
		All_families[x].append(File_Med_AllFams_Dct[y][name])
	Final_Dct[name] = All_families[x]
Final_Dct_Header = []
for x in samples:
	Final_Dct_Header.append(x)
Final_Dct['aaa_Family'] = Final_Dct_Header

### Section 9.) write output to file
outfile = open(MedOutName, 'w')
for key, value in sorted(Final_Dct.items()):
	outfile.write( str(key) + ',' + str(value).strip('[ | ]') + '\n' )

ifile.close()
mfile.close()		
