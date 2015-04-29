"""

OVERVIEW: 

Script to convert raw 16S sequencing data into OTU tables.  Acts on a directory containing a summary file and the raw data.
Outputs a directory with processing results.

"""

import preprocessing_16S as OTU
from optparse import OptionParser
import numpy as np
import os, sys
import os.path
from SummaryParser import *
import math
from string import ascii_lowercase
import multiprocessing as mp
import Formatting as frmt

# Read in arguments for the script

usage = "%prog -i INPUT_DIR -o OUTPUT_DIR_FULLPATH"
parser = OptionParser(usage)
parser.add_option("-i", "--input_dir", type="string", dest="input_dir")
parser.add_option("-o", "--output_dir", type="string", dest="output_dir")
parser.add_option("-p", "--primers_removed", dest="primers_removed", default='False')
parser.add_option("-b", "--split_by_barcodes", dest="split_by_barcodes", default='False')
(options, args) = parser.parse_args()

if( not options.input_dir ):
    parser.error("No data directory specified.")

# Parse summary file
summary_file = options.input_dir + '/summary_file.txt'
summary_obj = SummaryParser(summary_file)
summary_obj.ReadSummaryFile()
dataset_ID = summary_obj.datasetID

# If no output directory specified, default to $home/proc/
homedir = os.getenv("HOME")
if( not options.output_dir ):
    print "No output directory name specified.  Writing to " + homedir + "/proc/ by default."
    options.output_dir = homedir + '/proc/' + dataset_ID + '_proc_16S'

# Make a directory for the 16S processing results 
working_directory = options.output_dir
os.system('mkdir ' + working_directory)

# Extract file locations
raw_fastq_data_file = options.input_dir + '/' + summary_obj.attribute_value_16S['RAW_FASTQ_FILE']
primers_file = options.input_dir + '/' + summary_obj.attribute_value_16S['PRIMERS_FILE']
barcodes_map = options.input_dir + '/' + summary_obj.attribute_value_16S['BARCODES_MAP']

# Pipe stdout and stderr to logfiles in the new directory
sys.stdout = open(working_directory + '/stdout_proc_16S.log','w')
sys.stderr = open(working_directory + '/stderr_proc_16S.log','w')

# Construct output filenames from dataset ID
fastq_trimmed_qual = working_directory + '/' + dataset_ID + '.raw_trimmed_qual.fastq'
fasta_trimmed = working_directory + '/' + dataset_ID + '.raw_trimmed.fasta'
fastq_trimmed_length = working_directory + '/' + dataset_ID + '.raw_length_trimmed.fastq'
fastq_trimmed_primers = working_directory + '/' + dataset_ID + '.raw_primers_trimmed.fastq'
fastq_split_by_barcodes = working_directory + '/' + dataset_ID + '.raw_split_by_barcodes.fastq'
fasta_dereplicated = working_directory + '/' + dataset_ID + '.raw_dereplicated.fasta'
OTU_sequences_fasta = working_directory + '/' + dataset_ID + '.otu_seqs.fasta'
OTU_sequences_table = working_directory + '/' + dataset_ID + '.otu_seqs.table'
OTU_clustering_results = working_directory + '/' + dataset_ID + '.otu_clustering.tab'
OTU_database = working_directory + '/' + dataset_ID + '.otu_database'
OTU_table = working_directory + '/' + dataset_ID + '.otu_table'

# Get ASCII encoding of FASTQ files
try:
    encoding = summary_obj.attribute_value_16S['ASCII_ENCODING']
except:
    encoding = ''

if(encoding == "ASCII_BASE_33"):
    print "ASCII 33 encoding for quality scores specified."
    ascii_encoding = 33
elif(encoding == "ASCII_BASE_64"):
    print "ASCII 64 encoding for quality scores specified."
    ascii_encoding = 64
else:
    print "No ASCII encoding specified in the summary file for the quality scores in the FASTQ file.  Using ASCII 64 as default."
    ascii_encoding = 64

# Parallel steps:
#       1. split fastq into chunks
#       2. sort by barcodes, remove primers, and trim, and convert to fasta format
#       3. recombine into a single fasta file before dereplicating

# Step 1.1 - split file into 200000 line chunks
os.chdir(working_directory)
os.system('split -l 200000 ' + raw_fastq_data_file)

# Step 1.2 - get split filenames
split_filenames = []
for c1 in ascii_lowercase:
    for c2 in ascii_lowercase:
        filename = 'x'+c1+c2
        if(os.path.isfile(filename)):
            split_filenames.append(filename)
if len(split_filenames) == 0:
    split_filenames = [raw_fastq_data_file]

# Check whether samples need to be split by barcodes and primers need to be removed
#if (options.split_by_barcodes == 'True' and options.primers_removed == 'True'):
#    # Copy the raw file into processed folder and call it trimmed by primers
#    cmd_str = 'cp raw_fastq_data_file ' + fastq_trimmed_primers
#    os.system(cmd_str)

# Step 2 - loop through these split files and launch parallel threads
output = mp.Queue()

# Step 2.1 - split by barcodes
if (options.split_by_barcodes == 'False'):
    mode = summary_obj.attribute_value_16S['BARCODES_MODE']
    barcode_processes = [mp.Process(target=OTU.split_by_barcodes, args=(split_filename, split_filename + '.sb', barcodes_map, mode)) for split_filename in split_filenames]
    for p in barcode_processes:
        p.start()
    for p in barcode_processes:
        p.join()
    split_filenames = [f + '.sb' for f in split_filenames] 


# Step 2.2 - remove primers
if (options.primers_removed == 'False'):
    # Remove primers
    primer_processes = [mp.Process(target=OTU.remove_primers, args=(split_filename, split_filename + '.pt', primers_file)) for split_filename in split_filenames]
    for p in primer_processes:
        p.start()
    for p in primer_processes:
        p.join()
    split_filenames = [f + '.pt' for f in split_filenames] 



# Step 2.3 - trim with quality filter
quality_trim_processes = [mp.Process(target=OTU.trim_quality, args=(split_filename, split_filename + '.qt', ascii_encoding)) for split_filename in split_filenames]
for p in quality_trim_processes:
    p.start()
for p in quality_trim_processes:
    p.join()
split_filenames = [f + '.qt' for f in split_filenames] 


# Step 2.4 - trim to uniform length of 101
length = 101
length_trim_processes = [mp.Process(target=OTU.trim_length, args=(split_filename, split_filename + '.lt', length, ascii_encoding)) for split_filename in split_filenames]
for p in length_trim_processes:
    p.start()
for p in length_trim_processes:
    p.join()
split_filenames = [f + '.lt' for f in split_filenames] 

# Step 2.5 - convert to FASTA format
fasta_conversion_processes = [mp.Process(target=frmt.fastq2fasta, args=(split_filename, split_filename + '.fasta')) for split_filename in split_filenames]
for p in fasta_conversion_processes:
    p.start()
for p in fasta_conversion_processes:
    p.join()
split_filenames = [f + '.fasta' for f in split_filenames] 


# Step 3 - Recombine into a single fasta file
if len(split_filenames)>1:
    cat_str = ['cat']
    for filename in split_filenames:
        cat_str.append(filename)
    cat_str = ' '.join(cat_str)
    cat_str = cat_str + ' > ' + fasta_trimmed
    #cat_str = cat_str + ' > ' + fastq_trimmed_length
    
    # Recombine
    os.system(cat_str)


# Dereplicate and sort by size
try:
    separator = summary_obj.attribute_value_16S['BARCODES_SEPARATOR']
except:
    separator = '_'
OTU.dereplicate_and_sort(fasta_trimmed, fasta_dereplicated, OTU_database, separator)


# Remove chimeras
OTU.remove_chimeras_and_cluster_OTUs(fasta_dereplicated, OTU_sequences_fasta, OTU_sequences_table, OTU_clustering_results)


################
#
# Obtain Greengenes reference IDs for dereplicated sequences
#
################

alignment_results = working_directory + '/gg_alignments.aln'
cmd_str = 'usearch8 -usearch_local ' + OTU_sequences_fasta + ' -db /home/ubuntu/gg_13_5_otus/rep_set/97_otus.fasta -strand both -id 0.97 -alnout ' + alignment_results
os.system(cmd_str)


# Extract alignment dictionary
OTU_GG_dict = OTU.parse_alignment(alignment_results)

# Build 2 OTU tables - one complete, one only with OTUs that matched a Greengene sequence (for use in PiCRUST) 
OTU.build_OTU_table(OTU_sequences_fasta, OTU_database, OTU_table, OTU_GG_dict)

# Convert OTU tables to classic format
OTU_table_classic = OTU_table + '.classic'
OTU_table_gg = OTU_table + '.gg'
OTU_table_gg_classic = OTU_table_gg + '.classic'
frmt.convert_OTU_to_classic_dense_format(OTU_table, OTU_table_classic)
frmt.convert_OTU_to_classic_dense_format(OTU_table_gg, OTU_table_gg_classic)

# Convert to BIOM format
#OTU_table_biom = working_directory + '/' + dataset_ID + '.biom'
#OTU_table_biom_gg = working_directory + '/' + dataset_ID + '.gg.biom'
#frmt.build_OTU_table_biom(OTU_table + '.classic', OTU_table_biom, dataset_ID)
#frmt.build_OTU_table_biom(OTU_table + '.gg.classic', OTU_table_biom_gg, dataset_ID)

# Transfer to QIIME server for PiCRUST
cmd_str = 'scp -i /home/ubuntu/kys/ec2_private_key.pem ' + OTU_table_gg_classic + ' ubuntu@52.4.197.174:/home/ubuntu/picrust_inbox/.'
os.system(cmd_str)



###################
#
#  Final check
#
###################


# Check for file size greater than zero - add more thorough check eventually
otu_proc_success = False
if(os.stat(OTU_table).st_size > 0 and os.stat(OTU_sequences_fasta).st_size > 0 and os.stat(OTU_sequences_table).st_size > 0):
    otu_proc_success = True

# Processing complete - if successful, update summary file and write.  Otherwise, leave untouched and exit.
if(otu_proc_success == True):
    summary_obj.attribute_value_16S['OTU_TABLE'] = OTU_table
    summary_obj.attribute_value_16S['OTU_SEQUENCES_FASTA'] = OTU_sequences_fasta
    summary_obj.attribute_value_16S['OTU_SEQUENCES_TABLE'] = OTU_sequences_table
    summary_obj.attribute_value_16S['PROCESSED'] = 'True'
    summary_obj.WriteSummaryFile()
    print "Successfully processed 16S data!  Summary file has been updated."
else:
    print "Failed to process 16S data."


