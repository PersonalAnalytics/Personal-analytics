"""

OVERVIEW: 

Module to convert various file formats on the platform.

"""

import os, sys
import util
import numpy as np


def fasta2table(fastaIn, tableOut):
    # Converts a set of fasta sequences into a table format with the first column
    # corresponding to label lines beginning with > and the second column to the sequence.
    keep = {}
    seqs = {}

    for [otu_number, seq] in util.iter_fst(fastaIn):
        otu_number = otu_number[1:]
        keep[otu_number] = 1
        seqs[otu_number] = seq

    # Sort and organize into a new tab-delimited file with OTU_ID and Sequence as columns
    fid = open(tableOut,'w')
    headerline = "OTU_ID" + '\t' + 'Sequence'
    fid.write(headerline+'\n')
    for otu_number in keep:
        line = str(otu_number) + '\t' + seqs[otu_number]
        fid.write(line+'\n')
    fid.close()
    return None


def convert_OTU_table_to_classic_dense_format(otu_table_file, output_otu_table_file):
    # Converts an OTU table output from the 16S preprocessing steps into a
    # biom-format compatible classic format of the form:
    # OTU_ID  sample1   sample2
    # OTU0    0         0
    # OTU1    34        0
    # 
    # Note:
    #      -OTU labels are in the form 1,2,3,4 and get reassigned to OTU1, OTU2, OTU3, etc.
    #
    with open(otu_table_file,'r') as fidin:
        otu_table_data = fidin.readlines()
        firstrow = otu_table_data[0].split('\t')
        OTU_labels = firstrow[1:]
        sample_labels = [otu_table_data[i].split('\t')[0] for i in range(1,len(otu_table_data))]
        nOTUs = len(OTU_labels)
        nSamples = len(sample_labels)
        # Load OTU table row major order
        OTU_table_old = np.zeros((nSamples, nOTUs))
        for i in range(1,nSamples+1):
            OTU_table_old[i-1,:] = otu_table_data[i].split('\t')[1:]
        # Write transposed OTU table row major order
        OTU_table_t = np.transpose(OTU_table_old)
        with open(output_otu_table_file,'w') as fidout:
            fidout.write("OTU_ID" + '\t' + '\t'.join(sample_labels) + '\n')
            for i in range(nOTUs):
                tmpline = [str(int(OTU_table_t[i,j])) for j in range(len(OTU_table_t[i,:]))]
                line = OTU_labels[i] + '\t' + '\t'.join(tmpline)
                fidout.write(line + '\n')
