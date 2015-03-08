# usage: python upload_to_S3.py -f filename -b bucket_ARN -k filekey

"""

OVERVIEW:

This script uploads a given input file to a specified S3 bucket and give it a specified file key.  If no file key is specified, it will just use the same filename.

"""

# Written by Thomas Gurry 
# Updated: 02/22/2015

import boto
import boto.s3
import sys
from boto.s3.key import Key
from optparse import OptionParser

# Read in arguments for the script

usage = "%prog -f INPUT_FILE -b DESTINATION_BUCKET_ARN -k FILEKEY"
parser = OptionParser(usage)
parser.add_option("-f", "--input_file", type="string", dest="inputfile")
parser.add_option("-b", "--bucket_arn", type="string", dest="bucket_arn")
parser.add_option("-k", "--file_key", type="string", dest="fileKey")
(options, args) = parser.parse_args()

if( not options.inputfile ):
    parser.error("No file specified.")
if( not options.bucket_arn ):
    parser.error("No bucket Amazon Resource Name (ARN) specified.")
if( not options.fileKey ):
    fileKey = options.inputfile.split('/')
    fileKey = fileKey[len(fileKey)-1]
else:
    fileKey = options.fileKey


# Connect to S3
AWS_ACCESS_KEY_ID = XXXXXXXXXX
AWS_SECRET_ACCESS_KEY = XXXXXXXXXX
bucket_name = options.bucket_arn
conn = boto.connect_s3(AWS_ACCESS_KEY_ID,AWS_SECRET_ACCESS_KEY)
bucket = conn.create_bucket(bucket_name,location=boto.s3.connection.Location.DEFAULT)

print 'Uploading %s to Amazon S3 bucket %s' % \
    (options.inputfile, options.bucket_arn)

def percent_cb(complete,total):
    sys.stdout.write('.')
    sys.stdout.flush()

k = Key(bucket)
k.key = fileKey
k.set_contents_from_filename(options.inputfile,
                             cb=percent_cb, num_cb=10)

