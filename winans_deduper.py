#!/usr/bin/env python

import re
import os
import argparse
from subprocess import call

def get_args():
        parser = argparse.ArgumentParser(prog='deduper', description="Remove PCR duplicates from SAM file", add_help=False)
        parser.add_argument('-f', '--file', help = "input SAM file", required = True)
        parser.add_argument('-u', '--umi', help = "file containing UMIs", required = True)
        parser.add_argument('-p', '--paired', action = 'store_true', help = "paired-end data")
        return parser.parse_args() 
args = get_args()
in_sam = args.file
umi_file = args.umi
paired = args.paired

#extract input filename for output files
filename = os.path.basename(in_sam).split('.')[0]

if paired is True:
    print("Error: this program not compatible with paired-end data at this time")
    exit()

def position_fw(line):
    """Takes line of SAM file mapped to forward strand and, if soft clipped, returns soft-clipping adjusted postion, otherwise returns POS"""
    pos = int(line.split('\t')[3])
    CIGAR = line.split('\t')[5]
    if "S" in CIGAR:
        left_sc = re.match('^\d+S', CIGAR)
        if left_sc is not None:
            new_pos = pos - int(left_sc.group().split('S')[0])
            return(new_pos)
        else:
            return(pos)    
    else: 
        return(pos)
    
def position_rv(line):
    """Takes line of SAM file mapped to reverse strand and adds matches, mismatches, deletions, introns, and right soft-clipping to POS, returns adjusted position"""
    pos = int(line.split('\t')[3])
    CIGAR = line.split('\t')[5]
    adj = [int(i[:-1]) for i in re.findall("\d+S$|\d+M|\d+D|\d+N", CIGAR)]
    new_pos = pos + sum(adj)
    return(new_pos)


#sort input sam file by UMI, GNAME, POS w/bash script
call(['bash', 'sort.sh', in_sam])

sorted_sam = 'sorted_temp.sam'

#initialize dictionary with values = UMIs and keys = counts
with open(umi_file, 'rt') as umis:
    umi_dict = {}
    for umi in umis:
        umi_dict[umi.strip()] = 0

#open output files
#misindexed for unmapped reads and those with unknown or low quality UMIs
deduped = open('%s_deduped.sam' % filename, 'w')
duplicate = open('%s_duplicates.sam' % filename, 'w')
misindexed = open('%s_misindexed.sam' % filename, 'w') 

#initialize counters
deduped_count = 0
duplicate_count = 0
qual_count = 0
read_count = 0

seq_info_last = []
with open(sorted_sam, 'rt') as sorted_sam:
    for i,line in enumerate(sorted_sam):
        if line.startswith('@'): #write out header lines
            deduped.write(line)
        else:
            read_count += 1
            flag = int(line.split('\t')[1])
            if ((flag & 4) == 4): #write unmapped reads to misindexed
                misindexed.write(line)
                qual_count += 1
            else:
                umi = re.search('[ATCGN]+$', line.split('\t')[0]).group()
                if umi not in umi_dict: #write reads w/unknown umis to misindexed
                    misindexed.write(line)
                    qual_count += 1
                else:
                    umi_dict[umi] += 1
                    rname = int(line.split('\t')[2])
                    if ((flag & 16) == 16): #mapped to reverse strand
                        pos = position_rv(line)
                        strand = "rv"
                    else:
                        pos = position_fw(line)
                        strand = "fw"
                    seq_info_current = [umi, rname, pos, strand]
                    if seq_info_current == seq_info_last:
                        duplicate.write(line)
                        duplicate_count += 1
                    else:
                        deduped.write(line)
                        deduped_count += 1
                    seq_info_last = seq_info_current[:]

#close output files
deduped.close()
duplicate.close()
misindexed.close()

#remove temp file
os.remove('sorted_temp.sam')

#write summary stats to file
with open('deduper_summary.txt', 'w') as fh:
    fh.write("-------------------------\n")
    fh.write("| DEDUPLICATION SUMMARY |\n")
    fh.write("-------------------------\n\n")
    fh.write("Number of unique reads: " + str(deduped_count) + '\n')
    fh.write("Percent unique reads: " + str(round(float(deduped_count)/float(read_count)*100,2)) + '%\n')
    fh.write("Number of duplicate reads: " + str(duplicate_count) + '\n')
    fh.write("Percent duplicate reads: " + str(round(float(duplicate_count)/float(read_count)*100,2)) + '%\n')
    fh.write("Number of misindexed reads: " + str(qual_count) + '\n')
    fh.write("Percent of misindexed reads: " + str(round(float(qual_count)/float(read_count)*100,2)) + '%\n')
    fh.write("Total reads: " + str(read_count) + '\n')
    fh.write("\n\nRead counts per UMI:\n")
    fh.write("--------------------\n")
    for umi, count in sorted(umi_dict.items(), key=lambda x: x[1], reverse = True):
        fh.write(str(umi) + '\t' + str(count) + '\n')
