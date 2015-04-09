"""

OVERVIEW: 

Python module for creating, updating, and interpreting individual dataset summary files and objects.  These summary files are tab-delimited, standard format files that are intended to standardize and facilitate interactions between the user, the data, and different computing requests associated with the data.  They contain they information about the dataset in addition to the location of dataset-specific files.  

Data from the summary file can be loaded into a Python object, modified, and the changes written back to file.  In this manner, information about the dataset can be efficiently accessed and manipulated.


A summary file object is created in Python as follows:

obj = SummaryParser('summary_file.txt')


An example summary file is found between the dotted lines:

=============================================
DATASET_ID	obio_testdata

#16S
DATASET_ID	obio_testdata
RAW_FASTQ_FILE	obio.raw.fastq
PRIMERS_FILE	obio.primers.lst
BARCODES_MAP	obio.barcodes.lst
BARCODES_MODE   1
METADATA_FILE	obio.meta.txt
PROCESSED	True
ASCII_ENCODING  ASCII_64
AMPLICON	V1-V2
KEYWORDS	fmt

#CYTOKINES

#TCR

#METABOLOMICS

=============================================

"""

import numpy as np
import os
import os.path
import sys

class SummaryParser():
    def __init__(self, summary_file):
        self.summary_file = summary_file
        self.datasetID = None

        # Initialize 16S attributes
        self.attribute_value_16S = {'PROCESSED': "N/A"}

        
    def Extract16SLines(self):
        # Description:  Extracts the line numbers pertaining to 16S
        with open(self.summary_file,'r') as summary_fid:
            all_lines = summary_fid.readlines()
            for i in range(len(all_lines)):
                line = all_lines[i].split()
                if(len(line)>0):
                    if(line[0] == "#16S_start"):
                        startline = i+1
                    if(line[0] == "#16S_end"):
                        endline = i
                        break
        return [startline, endline]
        

    def ReadSummaryFile(self):
        # Description:  Loads the summary file specified in the object.
        with open(self.summary_file,'r') as summary_fid:
            summary_file_lines = summary_fid.readlines()
            # Read dataset ID
            self.datasetID = summary_file_lines[0].split('\t')[1].rstrip('\n')

            # Read 16S attributes
            [startline, endline] = self.Extract16SLines()
            for i in range(startline, endline):
                line = summary_file_lines[i].split('\t')
                attribute = line[0]
                value = line[1].rstrip('\n')
                self.attribute_value_16S[attribute] = value

 
    def WriteSummaryFile(self):
        # Description :  Writes a tab-delimited summary file using the values currently in the summary object.    
        with open(self.summary_file,'w') as summary_fid:

            # Populate summary file
            summary_fid.write('DATASET_ID' + '\t' + self.datasetID + '\n')
            summary_fid.write('\n')

            # 16S portion
            summary_fid.write('#16S_start' + '\n')
            for [key, value] in self.attribute_value_16S.items():
                summary_fid.write(key + '\t' + value + '\n')
            summary_fid.write('#16S_end' + '\n')

