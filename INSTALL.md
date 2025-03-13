
# PeptideMiner
PeptideMiner neuropeptide mining pipeline



Readme.txt file for PeptideMiner


# 1. Installation
1.1. Create conda enviroment for PeptideMiner:

conda create -n PeptideMiner -c bioconda -c conda-forge -c defaults hmmer blast sqlite python=3.10
conda activate PeptideMiner

1.2 Installing fast36
https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36/
wget https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8i.tar.gz
tar xzf fasta-36.3.8i.tar.gz
cd fasta-36.3.8i/src
make -f ../make/Makefile.linux64

1.2 Installing signalp
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=4.1&packageversion=4.1g&platform=Linux
tar xzf signalp-4.1g.Linux.tar.gz
cd signalp-4.1
vi signalp 
 edit line 13 - path to installed signalp-4.1
 
git clone https://github.com/fteufel/signalp-6.0
pip install signalp-6.0/

