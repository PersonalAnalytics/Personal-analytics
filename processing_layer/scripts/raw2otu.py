"""

OVERVIEW: 

Script to convert raw 16S sequencing data into OTU tables.  Acts on a directory containing a summary file and the raw data.
Outputs a directory with processing results.

"""

from __future__ import print_function
from optparse import OptionParser
import numpy as np
import os, sys
import os.path
import math
from string import ascii_lowercase
import multiprocessing as mp
import ntpath
import preprocessing_16S as OTU
import Formatting as frmt
from CommLink import *
from SummaryParser import *
from Features import *
import pickle

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

# Pipe stdout and stderr to logfiles in the new directory
sys.stdout = open('/home/ubuntu/logs/stdout_' + dataset_ID + '_proc_16S.log','w')
sys.stderr = open('/home/ubuntu/logs/stderr_' + dataset_ID + '_proc_16S.log','w')
def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)

# Check for presence of metadata map - report if metadata is missing.
try:
    metadata_file = summary_obj.attribute_value_16S['METADATA_FILE']
except:
    metadata_file = None
    warning("No metadata file found!!  This will cause problems downstream...")

# If no output directory specified, default to $home/proc/
homedir = os.getenv("HOME")
if( not options.output_dir ):
    print("No output directory name specified.  Writing to " + homedir + "/proc/ by default.")
    options.output_dir = homedir + '/proc/' + dataset_ID + '_proc_16S'

# Make a directory for the 16S processing results 
working_directory = options.output_dir
try:
    os.system('mkdir ' + working_directory)
except:
    print("Processing directory for this dataset already exists.  Overwriting its contents.")

# Extract file locations
try:
    raw_data_file = options.input_dir + '/' + summary_obj.attribute_value_16S['RAW_FASTQ_FILE']
    raw_file_type = 'FASTQ'
except:
    print("No raw FASTQ file found.  Checking for raw FASTA.")
    try:
        raw_data_file = options.input_dir + '/' + summary_obj.attribute_value_16S['RAW_FASTA_FILE']
        raw_file_type = 'FASTA'
    except:
        print("No raw FASTA file found either.  Check your raw data.")

primers_file = options.input_dir + '/' + summary_obj.attribute_value_16S['PRIMERS_FILE']
barcodes_map = options.input_dir + '/' + summary_obj.attribute_value_16S['BARCODES_MAP']

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
    print("ASCII 33 encoding for quality scores specified.")
    ascii_encoding = 33
elif(encoding == "ASCII_BASE_64"):
    print ("ASCII 64 encoding for quality scores specified.")
    ascii_encoding = 64
else:
    print ("No ASCII encoding specified in the summary file for the quality scores in the FASTQ file.  Using ASCII 64 as default.")
    warning("No ASCII encoding specified in the summary file for the quality scores in the FASTQ file.  Using ASCII 64 as default.")
    ascii_encoding = 64


# Parallel steps:
#       1. split fastq into chunks
#       2. demultiplex (sort by barcodes), remove primers, and trim, and convert to fasta format
#       3. recombine into a single fasta file before dereplicating

# Step 1.1 - get raw data filesize, then split into ~10Mb pieces (100000 lines) if smaller than 100 Mb, or into ~100Mb pieces (1000000 lines) otherwise.  Can optimize this eventually
# to split according to the number of cpus.
rawfilesize = os.path.getsize(raw_data_file)

# Step 1.1 - split file into 1000000 line (~100Mb) chunks
os.chdir(working_directory)
if(rawfilesize < 2e8):
    os.system('split -l 100000 ' + raw_data_file)
else:
    os.system('split -l 1000000 ' + raw_data_file)

# Step 1.2 - get split filenames
split_filenames = []
for c1 in ascii_lowercase:
    for c2 in ascii_lowercase:
        filename = 'x'+c1+c2
        if(os.path.isfile(filename)):
            split_filenames.append(filename)
if len(split_filenames) == 0:
    split_filenames = [raw_data_file]

# Check whether samples need to be split by barcodes and primers need to be removed
if (options.split_by_barcodes == 'True' and options.primers_removed == 'True'):
    # Copy the raw file into processed folder and call it trimmed by primers
    cmd_str = 'cp ' + raw_data_file + ' ' + fastq_trimmed_primers
    os.system(cmd_str)

# Step 2 - loop through these split files and launch parallel threads as a function of the number of CPUs
cpu_count = mp.cpu_count()

# Step 2.1 - demultiplex, i.e. sort by barcode
if (options.split_by_barcodes == 'False'):
    mode = summary_obj.attribute_value_16S['BARCODES_MODE']
    pool = mp.Pool(cpu_count)
    filenames = split_filenames
    newfilenames = [f + '.sb' for f in filenames]
    barcodes_map_vect = [barcodes_map]*len(filenames)
    mode_vect = [mode]*len(filenames)
    pool.map(OTU.split_by_barcodes, zip(filenames, newfilenames, barcodes_map_vect, mode_vect))
    pool.close()
    pool.join()
    split_filenames = [f + '.sb' for f in split_filenames] 


# Step 2.2 - remove primers
if (options.primers_removed == 'False'):
    pool = mp.Pool(cpu_count)
    filenames = split_filenames
    newfilenames = [f + '.pt' for f in filenames]
    primers_vect = [primers_file]*len(filenames)
    pool.map(OTU.remove_primers, zip(filenames, newfilenames, primers_vect))
    pool.close()
    pool.join()
    split_filenames = [f + '.pt' for f in split_filenames] 


# Step 2.3 - trim with quality filter
if (raw_file_type == "FASTQ"):
    pool = mp.Pool(cpu_count)
    filenames = split_filenames
    newfilenames = [f + '.qt' for f in filenames]
    ascii_vect = [ascii_encoding]*len(filenames)
    pool.map(OTU.trim_quality, zip(filenames, newfilenames, ascii_vect))
    pool.close()
    pool.join()
    split_filenames = [f + '.qt' for f in split_filenames] 


# Step 2.4 - trim to uniform length of 101
length = 101
pool = mp.Pool(cpu_count)
filenames = split_filenames
newfilenames = [f + '.lt' for f in filenames]
length_vect = [length]*len(filenames)
ascii_vect = [ascii_encoding]*len(filenames)
    
if (raw_file_type == "FASTQ"):
    pool.map(OTU.trim_length_fastq, zip(filenames, newfilenames, length_vect, ascii_vect))
    pool.close()
    pool.join()
    split_filenames = [f + '.lt' for f in split_filenames] 
else:
    pool.map(OTU.trim_length_fasta, zip(filenames, newfilenames, length_vect))
    pool.close()
    pool.join()
    split_filenames = [f + '.lt' for f in split_filenames] 


# Step 2.5 - convert to FASTA format
if (raw_file_type == "FASTQ"):
    pool = mp.Pool(cpu_count)
    filenames = split_filenames
    newfilenames = [f + '.fasta' for f in filenames]
    pool.map(frmt.fastq2fasta, zip(filenames, newfilenames))
    pool.close()
    pool.join()
    split_filenames = [f + '.fasta' for f in split_filenames] 

# Step 2.6 - renumber sequences IDs to be consistent across files
try:
    separator = summary_obj.attribute_value_16S['BARCODES_SEPARATOR']
except:
    separator = '_'
OTU.renumber_sequences(split_filenames, separator)

# Step 3 - Recombine into a single fasta file
if len(split_filenames)>1:
    cat_str = ['cat']
    for filename in split_filenames:
        cat_str.append(filename)
    cat_str = ' '.join(cat_str)
    cat_str = cat_str + ' > ' + fasta_trimmed    
    # Recombine
    os.system(cat_str)
else:
    os.system('cp ' + split_filenames[0] + ' ' + fasta_trimmed)

# Dereplicate sequences into a list of uniques for clustering
OTU.dereplicate_and_sort(fasta_trimmed, fasta_dereplicated, OTU_database, '_')

# Remove chimeras and cluster OTUs
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


#################################################
#
# Put all results and metadata file in a single folder
#
#################################################

dataset_folder = dataset_ID + '_results'
try:
    os.system('mkdir ' + dataset_folder)
except:
    print('Results directory already exists.  Overwriting its contents.')
os.system('cp ' + OTU_table_gg_classic + ' ' + dataset_folder + '/.')
os.system('cp ' + OTU_table_classic + ' ' + dataset_folder + '/.')
os.system('cp ' + OTU_sequences_table + ' ' + dataset_folder + '/.')
os.system('cp ' + OTU_sequences_fasta + ' ' + dataset_folder + '/.')
if metadata_file is not None:
    os.system('cp ' + options.input_dir + '/' + metadata_file + ' ' + dataset_folder + '/.')

# Put the summary file in the folder and change the summary file path to its new location
os.system('cp ' + summary_file + ' ' + dataset_folder + '/.')
summary_obj.summary_file = dataset_folder + '/summary_file.txt'
summary_obj.attribute_value_16S['OTU_TABLE_CLOSED_REF'] = ntpath.basename(OTU_table_gg_classic)
summary_obj.attribute_value_16S['OTU_TABLE_OPEN_REF'] = ntpath.basename(OTU_table_classic)
summary_obj.attribute_value_16S['OTU_SEQUENCES_TABLE'] = ntpath.basename(OTU_sequences_table)
summary_obj.attribute_value_16S['OTU_SEQUENCES_FASTA'] = ntpath.basename(OTU_sequences_fasta)
try:
    summary_obj.attribute_value_16S['METADATA_FILE'] = ntpath.basename(metadata_file)
except:
    summary_obj.attribute_value_16S['METADATA_FILE'] = "None"
summary_obj.attribute_value_16S['PROCESSED'] = 'True'
summary_obj.WriteSummaryFile()

# Transfer to PiCRUST server and wait for results
cl = CommLink('proc')
results_folder = dataset_ID + '_results'
test = cl.launch_proc_listener(dataset_folder, results_folder)

# Move results from inbox to results folder
processing_results_dir = '/home/ubuntu/processing_results'
os.system('mv /home/ubuntu/inbox/' + results_folder + ' ' + processing_results_dir + '/.')
os.chdir(os.path.join(processing_results_dir, results_folder))

# Extract features
features = Features('summary_file.txt')
features.LoadOTUtable()
features.LoadPredictedMetagenome()
metapredL1 = os.path.join('/home/ubuntu/processing_results', results_folder, 'picrust_results/CRC_Zhao_2012.L1.biom')
metapredL2 = os.path.join('/home/ubuntu/processing_results', results_folder, 'picrust_results/CRC_Zhao_2012.L2.biom')
metapredL3 = os.path.join('/home/ubuntu/processing_results', results_folder, 'picrust_results/CRC_Zhao_2012.L3.biom')
features.LoadPredictedMetagenome(metapredL1)
features.LoadPredictedMetagenome(metapredL2)
features.LoadPredictedMetagenome(metapredL3)
features.LoadPhylogeneticFeatures()

# Pickle features
with open('pickled_features.pkl', 'wb') as fid:
    pickle.dump(features, fid)


###################
#
#  Final check - implement a better check in the future
#
###################


# Check for file size greater than zero - add more thorough check eventually
otu_proc_success = False
if(os.stat(OTU_table).st_size > 0 and os.stat(OTU_sequences_fasta).st_size > 0 and os.stat(OTU_sequences_table).st_size > 0):
    otu_proc_success = True

# Processing complete - if successful, update summary file and write.  Otherwise, leave untouched and exit.
if(otu_proc_success == True):
    print("Successfully processed 16S data!  Summary file has been updated.")
else:
    print("Failed to process 16S data.")


