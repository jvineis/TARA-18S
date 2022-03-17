#!/usr/bin/env python

"""
    Read all fasta files and build a sorted OTU contingency
    table.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>, edited by Joe Vineis <jvineis@gmail.com>"
__date__ = "2022/03/15"
__version__ = "1.0"

import os
import re
import sys
import operator
import argparse

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#

parser = argparse.ArgumentParser(description='''This script takes the output returned from swarmv3.1 - when run like this "swarm -d 1 -f -t 40 -z pooled-samples-derep.fa -i pooled-samples-derep-struct.txt -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt"''')
parser.add_argument('-repfa', help = "The fasta of represenataive sequences for each swarm (-w flag) and then sorted using vsearch")
parser.add_argument('-stats' , help = "The swarm stats file (-s flag)")
parser.add_argument('-swarms', help = "The swarms table output of the -o flag in swarm")
parser.add_argument('-l' , help = "A single column list of the names of the dereplicated fasta files-should look something like this ERR562370-primer-derep - and the corresponding fasta file should look like this ERR562370-primer-derep.fa")
args = parser.parse_args()

def representatives_parse():
    """
    Get seed sequences.
    """
    separator = ";size="
    representatives_file = args.repfa
    representatives = dict()
    with open(representatives_file, "r") as representatives_file:
        for line in representatives_file:
            if line.startswith(">"):
                amplicon = line.strip(">;\n").split(separator)[0]
            else:
                representatives[amplicon] = line.strip()

    return representatives

def stats_parse():
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = args.stats
    stats = dict()
    seeds = dict()
    with open(stats_file, "r") as stats_file:
        for line in stats_file:
            cloud, mass, seed, seed_abundance = line.strip().split(separator)[0:4]
            stats[seed] = int(mass)
            seeds[seed] = (int(seed_abundance), int(cloud))
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.items(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds

def swarms_parse():
    """
    Map OTUs.
    """
    separator = "_[0-9000000]+|;size=[0-9000000]+;?| "  # parsing of abundance annotations
    swarms_file = args.swarms
    swarms = dict()
    with open(swarms_file, "r") as swarms_file:
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            seed = amplicons[0]
            swarms[seed] = [amplicons]

    return swarms


def fasta_parse():
    """
    Map amplicon ids, abundances and samples. These must be dereplicated fastas. 
    """
    separator = ";size="
    fasta_files = open(args.l, 'r')
    samples = dict()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        fa = fasta_file.strip()+".fa"
        sample = fasta_file.split("-")[0]
        samples[sample] = samples.get(sample, 0) + 1
        with open(fa, "r") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon, abundance = line.strip(">;\n").split(separator)
                    abundance = int(abundance)
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        # deal with duplicated samples
                        amplicons2samples[amplicon][sample] = amplicons2samples[amplicon].get(sample, 0) + abundance
    # deal with duplicated samples
    duplicates = [sample for sample in samples if samples[sample] > 1]
    if duplicates:
        print("Warning: some samples are duplicated", file=sys.stderr)
        print("\n".join(duplicates), file=sys.stderr)
    samples = sorted(samples.keys())
    
    return amplicons2samples, samples


def print_table(representatives, stats, sorted_stats,
                swarms, amplicons2samples,
                samples, seeds):
    """
    Export results.
    """
    # Print table header
    print("OTU", "total", "cloud",
          "amplicon", "length", "abundance",
          "spread",
          "sequence",
          "\t".join(samples),
          sep="\t", file=sys.stdout)

    # Print table content
    i = 1
    for seed, abundance in sorted_stats:
        sequence = representatives[seed]
        occurrences = dict([(sample, 0) for sample in samples])
        for amplicons in swarms[seed]:
            for amplicon in amplicons:
                for sample in samples:

                    occurrences[sample] += amplicons2samples[amplicon].get(sample, 0)
        spread = len([occurrences[sample] for sample in samples if occurrences[sample] > 0])
        sequence_abundance, cloud = seeds[seed]

        # output
        print(i, abundance, cloud,
              seed, len(sequence), sequence_abundance,
              spread, sequence,
              "\t".join([str(occurrences[sample]) for sample in samples]),
              sep="\t", file=sys.stdout)
        i += 1

    return

def main():
    """
    Read all fasta files and build a sorted OTU contingency table.
    """
    # Parse the representative sequences for each swarm
    representatives = representatives_parse()

    # Parse stats
    stats, sorted_stats, seeds = stats_parse()

    # Parse swarms
    swarms = swarms_parse()

    # Parse fasta files
    amplicons2samples, samples = fasta_parse()

    # Print table header
    print_table(representatives, stats, sorted_stats, swarms,
                amplicons2samples, samples,
                seeds)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
