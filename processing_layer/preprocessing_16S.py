"""

OVERVIEW: 

Python module for pre-processing FASTQ data into OTU tables.

"""

import numpy as np
import sys
import os, sys
import util

def length_stats_fastq(fastq_in):
    # Returns 5th percentile of read length for a 100000 sample from a FASTQ file.
    iter_seq = util.iter_fsq
    x = []
    counter = 0
    for record in iter_seq(fastq_in):
        sid, seq = record[:2]
        counter = counter + 1
        if(counter > 100000):
            break
        x.append(len(seq))
    x = np.array(x)
    return np.percentile(x,5)

def trim_quality(fastq_in, fastq_out):
    # Finds maximum Q-score cut-off where 95% of reads are over 200 bases, if possible.  If not, selects Q=5 and returns a warning.
    Qvals = range(5,11)
    bestQ = 0

    for Q in Qvals:
        # Trim to Q
        str1 = 'usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(Q) + ' -fastq_ascii 64 -fastqout ' + fastq_out
        os.system(str1)
        # Check length distribution
        fifthPercentile = length_stats_fastq(fastq_out)
        if(fifthPercentile >= 200):
            bestQ = Q
            bestFifthPercentile = fifthPercentile
        else:
            break

    if (bestQ == 0):
        print "[[ Quality trimming ]] ERROR!!  Could not obtain 95% of reads over 200 base pairs with quality score cut-off of at least 5.  Check sequencing quality of the data."
        return None
    else:
        str1 = 'usearch8 -fastq_filter ' + fastq_in + ' -fastq_truncqual ' + str(bestQ) + ' -fastq_ascii 64 -fastqout ' + fastq_out
        os.system(str1)
        print "[[ Quality trimming ]] Input file: " + fastq_in
        print "[[ Quality trimming ]] Trimmed sequences with quality cut-off of Q = " + str(bestQ) + "."  
        print "[[ Quality trimming ]] 5th percentile of read length: " + str(bestFifthPercentile) + "."
        print "[[ Quality trimming ]] Trimmed file: " + fastq_out
        print "[[ Quality trimming ]] Complete."
        return None

def trim_length(fastq_in, fastq_out, length):
    # Trims to uniform length and filters by maximum expected error
    str1 = 'usearch8 -fastq_filter ' + fastq_in + ' -fastq_trunclen ' + str(length) + ' -fastq_maxee 0.25 -fastq_ascii 64 -fastqout ' + fastq_out
    os.system(str1)
    # Check for how many reads were thrown out
    statinfoIN = os.stat(fastq_in)
    input_filesize = float(statinfoIN.st_size)
    statinfoOUT = os.stat(fastq_out)
    output_filesize = float(statinfoOUT.st_size)
    percent_thrown_out = 100*(output_filesize / input_filesize)
    print "[[ Length trimming ]] Input file: " + fastq_in
    print "[[ Length trimming ]] Trimmed sequences with length less than " + str(length) + " and maximum expected error of 0.25."
    print "[[ Length trimming ]] Threw out " + str(percent_thrown_out) + " % of reads."
    print "[[ Length trimming ]] Complete."
    return None
 
def remove_primers(fastq_in, fastq_out, primers_file):
    # Remove primers from FASTQ file
    print "[[ Primer trimming ]] ..."
    os.system('python ~/scripts/1.remove_primers.py -q ' + fastq_in + ' -l ' + primers_file + ' -d 1 -o ' + fastq_out)
    print "[[ Primer trimming ]] Complete."
    return None


def split_by_barcodes(fastq_in, fastq_out, barcodes_map):
    # Split by barcodes
    print "[[ Splitting by barcodes ]] ..."
    os.system('python ~/scripts/2.split_by_barcodes.py -q ' + fastq_in + ' -b ' + barcodes_map + ' -B tab -d 1 --mode 1 -o ' + fastq_out)
    return None

def dereplicate_and_sort(fastq_in, fasta_out, OTU_database):
    # Dereplicate and sort sequences by size
    print "[[ Dereplicating and sorting ]] ..."
    os.system('python ~/scripts/3.dereplicate.py -q ' + fastq_in + " -s ';' -o " + OTU_database + ' -d ' + fasta_out)
    print "[[ Dereplicating and sorting ]] Complete."
    return None


def remove_chimeras_and_cluster_OTUs(fasta_in, fasta_out):
    # Remove chimeric sequences and then cluster OTUs with default similarity of 97%.
    print "[[ Removing chimeras and clustering OTUs ]] ..."
    os.system('usearch8 -cluster_otus ' + fasta_in + ' -otus ' + fasta_out + ' -sizein')
    print "[[ Removing chimeras and clustering OTUs ]] Complete."
    return None


def build_OTU_table(fasta_in, dereplication_map, fasta_out):
    # Builds an OTU table from a list of OTUs and a dereplication map
    print "[[ Building OTU table ]] ..."
    os.system('python ~/scripts/4.derep2counts.py --fst ' + fasta_in + ' --map ' + dereplication_map + ' --out ' + fasta_out)
    print "[[ Building OTU table ]] Complete."
    return None


fastq_raw = '/home/ubuntu/tmp_S3/obio.raw.fastq'
fastq_trimmed_qual = '/home/ubuntu/tmp_S3/obio.raw_trimmed_qual.fastq'
fastq_trimmed_length = '/home/ubuntu/tmp_S3/obio.raw_trimmed_length.fastq'
fastq_trimmed_primers = '/home/ubuntu/tmp_S3/obio.raw_trimmed_primers.fastq'
fastq_split_by_barcodes = '/home/ubuntu/tmp_S3/obio.fastq'
fasta_dereplicated = '/home/ubuntu/tmp_S3/obio_derep.fasta'
OTU_counts = '/home/ubuntu/tmp_S3/obio_otus.counts'
primers_file = '/home/ubuntu/tmp_S3/obio.primers.lst'
barcodes_map = '/home/ubuntu/tmp_S3/obio.barcodes.lst'
OTU_database = '/home/ubuntu/tmp_S3/obio.derep.map'

trim_quality(fastq_raw, fastq_trimmed_qual)
trim_length(fastq_trimmed_qual, fastq_trimmed_length, 260)
remove_primers(fastq_trimmed_length, fastq_trimmed_primers, primers_file) 
split_by_barcodes(fastq_trimmed_primers, fastq_split_by_barcodes, barcodes_map)
dereplicate_and_sort(fastq_split_by_barcodes, fasta_dereplicated, OTU_database)
remove_chimeras_and_cluster_OTUs(fasta_dereplicated, OTU_counts)

