#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''This script will combine the NODE-HITS.txt table from ASVs searched against the Euk 18S V9 Pr2 database two files will be created, one a count matrix, and the other a metadata file with the details of each SWARM''')
parser.add_argument('-n' , help = "The NODE-HITS.txt table")
parser.add_argument('-o' , help = "The prefix for the file names of the output files you like")
parser.add_argument('-r' , help = "The reference taxonomy file created from the reference fasta file") 
parser.add_argument('-s' , help = "The swarm count contingency table produced by ~/scripts/mu-swarms-to-ASVs-table-for-tarra.py")
parser.add_argument('-min', help = "The minimum number of sequences within a SWARM to be included in the output table")
args = parser.parse_args()

outfile_counts = open(args.o+'-counts.txt', 'w')
outfile_meta = open(args.o+'-metadata.txt','w')
node_hits = {}
for line in open(args.n,'r'):
    x = line.strip().split('\t')
    node_hits[x[0].split(";")[0]] = x[0:len(x)]

swarms = {}
for line in open(args.s, 'r'):
    x = line.strip().split('\t')
#    print(x[0:7])
    xasn = x[8:len(x)]
    if x[0] == "OTU":
        outfile_meta.write('\t'.join(x[0:7])+'\t')
        outfile_counts.write("OTU"+'\t'+'\t'.join(x[8:len(x)])+'\n')
    else:
        xl = []
        xasn = x[8:len(x)]
        for i in xasn:
    	    xl.append(int(i))
        if sum(xl) > int(args.min):
#            print(sum(xl), args.min)
            swarms[x[3]]=x[0:len(x)]


outfile_meta.write("taxid"+'\t'+"taxonomy"+'\t'+"tax1"+'\t'+"tax2"+'\t'+"tax3"+'\t'+"tax4"+'\t'+"tax5"+'\t'+"tax6"+'\t'+"tax7"+'\t'+"tax8"+'\t'+"tax9"+'\t'+"tax10"+'\n')
        
taxa_ref = {}
for line in open(args.r, 'r'):
    x = line.strip().split('\t')
    taxa_ref[x[0]] = x[0:len(x)]


for key in swarms.keys():
    if key in node_hits.keys():
        x = taxa_ref[node_hits[key][1]][1].split("|")
        outfile_meta.write("s_"+swarms[key][0]+'\t'+'\t'.join(swarms[key][1:7])+'\t'+'\t'.join(taxa_ref[node_hits[key][1]])+'\t'+'\t'.join(x)+'\n')
        outfile_counts.write("s_"+swarms[key][0]+'\t'+'\t'.join(swarms[key][8:len(swarms[key])])+'\n')
    else:
        outfile_meta.write("s_"+swarms[key][0]+'\t'+'\t'.join(swarms[key][1:7])+'\t'+"UNKNOWN"+'\t'+"UNKNOWN"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\t'+"unknown"+'\n')
        outfile_counts.write("s_"+swarms[key][0]+'\t'+'\t'.join(swarms[key][8:len(swarms[key])])+'\n')
                            
