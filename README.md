# PeptideMiner
PeptideMiner neuropeptide mining pipeline

Readme.txt file for PeptideMiner

PeptideMiner is a pipeline that uses profile hidden Markov models (profile-HMMs) of neuropeptide families to search query datasets for neuropeptide homologues. This document outlines the program requirements (section 1), how to do a test run (section 2) and the profile-HMMs it comes with (section 3). 


# 1. PeptideMiner requirements and installation
1) Programs that need to be installed for PeptideMiner to run include:

- signal-4.1

- python and conda (miniconda or anaconda)
    blast (makeblastdb and blastp), fasta36 and hmmer (hmmsearch) can be install via conda (bioconda) or installed as standalone

2) Create conda enviroment for PeptideMiner:

conda create -n PeptideMiner -c bioconda -c conda-forge -c defaults hmmer blast sqlite fasta3 configargparse numpy python=3.10
conda activate PeptideMiner

3) Installing signalp

https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=4.1&packageversion=4.1g&platform=Linux
tar xzf signalp-4.1g.Linux.tar.gz
cd signalp-4.1
vi signalp 
 edit line 13 - path to installed signalp-4.1

vi PeptideMiner.cfg
 edit line 12 - signalp_path : <path to installed signalp-4.1>/signalp

PeptideMiner is currently validated only for SignalP v4.1. 
Newer version of SignalP are available (https://services.healthtech.dtu.dk/) and while they should work with PeptideMiner, they have higher hardware and software requirements.   


# 2. Running PeptideMiner

1) The user creates a directory from where the program will run (working directory). All output files will be directed to the working directory

2) Copy the PeptideMiner.cfg file into this directory and update the PeptideMiner.cfg file with the user PATHS and search variables.

- User has to add path to signalp 4.1

3) Query databases must be in amino acids, fasta format and have the .fna extension. You can search multiple query databases in one search by putting them in the same directory.

4) Run the program from the command line from the working directory: 

    python PATH/TO/PeptideMiner/PeptideMiner.py --peptide_family <name of peptide family> --query <folder of query fasta files> 

  Parameters for command line :
    --peptide_family : Name of the peptide family. 
        Requires a folder in <data_dir>/01-known-seq with fasta file(s) of known sequences for this family     
        Requires at least one HMM profile (*.hmm) in <data_dir>/02-pHMM with <peptide_family> in its filename

    --query Folder to fasta files used for the HMMSearch query

  Parameters given in PeptideMiner.cfg or on command line (overwriting PeptideMiner.cfg):
    --cds-min_length : Minimum lenght of the protein-coding sequences CDS (Default 50)
    --signalp_cutoff : Cutoff of SignalP D-value [discrimination between singnal and non-signal peptides] (Default 0.41)
    --signalp_min_length : Minimum lenght of the signal peptides (Default 8)
    --signalp_path : Absolute path to signalp program
    --mature_min_length : Minimum lenght of mature peptids (Default 7)
    --mature_max_length : Maximum lenght of mature peptides (Default 15)
    --mature_evalue_cutoff : Cutoff of Fasta36 E-value alignment (Default 1)


# 3. Test PeptideMiner
1) Copy the Test.cfg from ./data/09-test/ to  chosen working directory

2) Adjust the signalp_path variable in the Test.cfg:

signalp_path : <path to signalp 4.1>/signalp

3) Run python PeptideMiner.py --config Test.cfg

4) Compare output in ./work/02-pipeline to output in ./data/09-test/02-pipeline


# 4. Available Data
1) Available profile-HMMs

PeptideMiner comes with five profile-HMMs :
insulin, natriuretic peptide, oxytocin, somatostatin, and tachykinin. 

profile-HMMs currently available:

Choose from:
	
  - test_OTVP
	
  - insulin
	
  - natriureticpeptides
	
  - oxytocin
	
  - tachykinin

It can search multiple query databases with more than one profile-HMM for the same neuropeptide family. To achieve this, PeptideMiner splits the file names of each profile-HMM in the data/02-pHMM directory at '_'. It matches the phrase preceding the first '_' to the  "neuropeptide_family" name provided in the config.txt and uses the profile-HMMs that match.

3.2. Using your own profile-HMM

When using your own profile-HMM the first word in the file name needs to match the "neuropeptide_family" specified in the config.txt file. The same applies to the directory where the known mature peptide sequences are stored in the data/01-known_seq directory. 

