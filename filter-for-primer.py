#!/usr/bin/env python

from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='''find the primer sequences for 18S used in the DeVargas paper in each sequence and create a new fasta containing only sequences that have both the fwd and reverse primer''')
parser.add_argument('--i', help='a merged fasta file')
parser.add_argument('--o', help='the name you want to give to your beautiful fasta that contains only the best sequences that contain the primer')
args=parser.parse_args()

outfile = open(args.o, 'w')

def filter_left_primers(records, Fprimer,Fprimer1,Rprimer):
    left_side_found_ids = []
    for record in records: # look for the fwd and rev primers in the sequence
        if record.seq.startswith(Fprimer) or record.seq.startswith(Fprimer1) or record.seq.startswith(Rprimer):
            left_side_found_ids.append(record)
    print("found this many good left", len(left_side_found_ids))
    return(left_side_found_ids)


def filter_right_primers(seq_list, Fprimer, Fprimer1, Rprimer):
    both_primers_found_ids = []
    for record in seq_list:
        if record.seq.endswith(Fprimer) or record.seq.endswith(Fprimer1) or record.seq.endswith(Rprimer):
            both_primers_found_ids.append(record)
    print("found this many right", len(both_primers_found_ids))
    return(both_primers_found_ids)
                                          
                               
original_reads = SeqIO.parse(args.i, "fasta")
left_trim = filter_left_primers(original_reads, "ttgtacacaccgccc","ccttccgcaggttcacctac","ccttctgcaggttcacctac")
both_found = filter_right_primers(left_trim, "gtaggtgaacctgcggaagg","gtaggtgaacctgcagaag", "gggcggtgtgtacaa")
SeqIO.write(both_found, args.o, "fasta")
