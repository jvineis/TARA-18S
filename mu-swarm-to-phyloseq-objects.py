#!/usr/bin/env python
from Bio import SeqIO
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description='''The script takes the output from SWARM and vsearch to create a taxonomy object and a matrix object that is compatible with phyloseq''')
parser.add_argument('-tax_ref', help = 'the reference taxonomy table that matches the ref fasta file searched agaist with a taxonomy string e.g. silva.tax')
parser.add_argument('-hits', help = 'the output from vsearch when run like this - vsearch --usearch_global reduced-node-fasta-min10-famsa.fa --db /groups/g454/blastdbs/gast_distributio\
ns/silva119.fa --blast6out NODE-HITS.txt --id 0.6  *** keep in mind that the --db is dependent on the gene of interest')
parser.add_argument('-med', help = 'the count or percent abundance matrix produced by mu-swarms-to-ASVs-table.py')
parser.add_argument('-fa', help = 'the fasta file of node hits produced by mu-swarms-to-ASVs-table.py - could be something like this "reduced-node-fasta-min10-famsa.fa"')
args = parser.parse_args()
 
def open_node_hits_to_dict(sample_name):
    sample_dict = {}
    with open(sample_name, 'rU') as f:
        for line in f:
            x = line.strip().split("\t")
            sample_dict[x[0]] = x[1]
    return sample_dict

def open_silva_database_to_dict(db_name):
    db_name_dict = {}
    with open(db_name, 'rU') as f:
        for line in f:
            x = line.rstrip().split('\t')
            db_name_dict[x[0]] = x[1:len(x)]
    return db_name_dict
 
def open_med_table_to_dict(med_table):
    med = {}
    with open(med_table, 'rU') as f:
        table_l = []
        for line in f:
            x = line.rstrip().split('\t')
            table_l.append(x)
        ## transpose the table
        t_table_l = [[table_l[j][i] for j in range(len(table_l))]for i in range(len(table_l[0]))]
        for line in t_table_l:
            med[line[0].split("|")[0]] = line[1:len(line)]
    return med

def open_NODE_REPRESENTATIVES(node_fasta):
    all_nodes = {}
    with open(node_fasta, 'rU') as f:
        x = SeqIO.parse(f, "fasta")
        for rec in x:
            all_nodes[rec.id] = rec.id
    return all_nodes
# Create each of the dictionaries using the defined functions above
 
node_hits_dict = open_node_hits_to_dict(args.hits)
tax_lookup_dict = open_silva_database_to_dict(args.tax_ref)
med_matrix_dict = open_med_table_to_dict(args.med)
all_nodes_dict = open_NODE_REPRESENTATIVES(args.fa)


#  Open and output file
output_tax = open('PHYLOSEQ-TAX-OBJECT.txt', 'w')
output_matrix = open('PHYLOSEQ-MATRIX-OBJECT.txt', 'w')

output_matrix.write("node"+'\t'+'\t'.join(med_matrix_dict['samples'])+'\n')

for key in med_matrix_dict.keys():
    print(key)

for key in all_nodes_dict.keys():
#    print all_nodes_dict[key].split("|")[0]
    if all_nodes_dict[key].split(";")[0] in med_matrix_dict.keys(): # added .split(";")[0] ####

        output_matrix.write(key.split(';')[0]+'\t'+'\t'.join(med_matrix_dict[all_nodes_dict[key].split("|")[0]])+'\n')

#Add header to output file
output_tax.write("node"+';'+"db_hit_id"+';'+"Kingdom"+";"+"Phylum"+';'+"Class"+';'+"Order"+';'+"Family"+';'+"Genus"+';'+"Species"+';'+"Environment"+";"+"ASV"+'\n')
 
# Add the node and taxonomy string to the output file
for key in node_hits_dict.keys():
    genome = node_hits_dict[key].split("|")[0]
    
    output_tax.write(str(key.split(";")[0])+';'+str(node_hits_dict[key].split(";")[0])+';'+';'.join(tax_lookup_dict[genome])+';'+str(key.split(";")[0].split("|")[0])+'\n')
    
#  Add the nodes that rec'd no hit in the  silva db to the taxonomic output.
for key in all_nodes_dict.keys():
    if key not in node_hits_dict.keys():
        output_tax.write(str(key.split(";")[0])+';'+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+str(key.split(";")[0].split("|")[0])+'\n')

output_tax.close()
output_matrix.close()
