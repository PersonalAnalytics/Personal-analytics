# usage: python download_from_S3.py -k fileKey -b bucket_ARN -f filename

"""

OVERVIEW:

This script downloads a given file key from a specified S3 bucket to a specified filename and path.

"""

# Written by Thomas Gurry 
# Updated: 02/22/2015

import boto
import sys, os
import boto.s3
from boto.s3.key import Key
from optparse import OptionParser

# Read in arguments for the script

usage = "%prog -k FILE_KEY -b BUCKET_ARN -f FILENAME"
parser = OptionParser(usage)
parser.add_option("-k", "--fileKey", type="string", dest="fileKey")
parser.add_option("-b", "--bucket_arn", type="string", dest="bucket_arn")
parser.add_option("-f", "--filename", type="string", dest="filename")
(options, args) = parser.parse_args()

if( not options.fileKey ):
    parser.error("No file key specified.")
if( not options.bucket_arn ):
    parser.error("No bucket Amazon Resource Name (ARN) specified.")
if( not options.filename ):
    parser.error("No filename or path specified as download location.")

fileKey = options.fileKey

# Connect to S3
AWS_ACCESS_KEY_ID = 'AKIAJ7BESVOEO2UKRRXQ'
AWS_SECRET_ACCESS_KEY = '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'
bucket_name = options.bucket_arn
conn = boto.connect_s3(AWS_ACCESS_KEY_ID,AWS_SECRET_ACCESS_KEY)
bucket = conn.create_bucket(bucket_name,location=boto.s3.connection.Location.DEFAULT)

# Download file

print 'Downloading %s from Amazon S3 bucket %s' % \
    (options.fileKey, options.bucket_arn)

def percent_cb(complete,total):
    sys.stdout.write('.')
    sys.stdout.flush()

k = Key(bucket)
k.key = fileKey

# Check if file exists locally, if not: download it
filename = options.filename
if not os.path.exists(filename):
    k.get_contents_to_filename(filename,cb=percent_cb, num_cb=10)

