#!/bin/bash
# Description: This script downloads and install all relevant packages to a naked ubuntu
# EC2 instance. Successful execution of this script will have created a machine that is
# equivalent to the latest processing layer EC2 AMI.
# Author: Thomas Gurry
#########################################
# General setup (text editors, etc. #
#########################################
mkdir ~/downloads
sudo chown ubuntu:ubuntu ~/downloads
sudo apt-get update
sudo apt-get install -y emacs24 --fix-missing
sudo apt-get upgrade gcc
sudo apt-get install gfortran


#################################
# Python setup #
#################################
# Get PIP
cd ~/downloads
wget https://bootstrap.pypa.io/get-pip.py
sudo chown ubuntu:ubuntu get-pip.py
python get-pip.py
sudo apt-get install -y build-essential python-dev python-setuptools \
python-scipy libatlas-dev libatlas3gf-base --fix-missing


# Make sure ATLAS is used to provide implementation of BLAS and LAPACK linear algebra routines
sudo update-alternatives --set libblas.so.3 \
/usr/lib/atlas-base/atlas/libblas.so.3
sudo update-alternatives --set liblapack.so.3 \
/usr/lib/atlas-base/atlas/liblapack.so.3

# Install numpy, scikit-learn, etc.
sudo pip install numpy==1.5.1  # Version required for PiCRUST
pip install --user --install-option="--prefix=" -U scikit-learn
pip install --user Biopython
pip install --user pandas
pip install --user boto
pip install --user bidict

#################################
# mySQL setup #
#################################
sudo apt-get install mysql-server mysql-client -y
sudo apt-get install python-mysqldb -y

###############
#  AWS setup  #
###############

sudo apt-get install awscli -y
sudo apt-get install s3cmd

##########################
#  Misc. bioinformatics  #
##########################
# sff2fastq
git clone git://github.com/indraniel/sff2fastq.git
cd sff2fastq
make
