import os
from libfunc import config,output

"""
BLAST module.
Support module for Step 7.
"""

def make(knownseqfile):
    #Create a ncbi db of knownseqfile
    blastdb = os.popen(f'{0}/makeblastdb -in {knownseqfile} -dbtype prot -hash_index -out {knownseqfile}')
    # blastdb = os.popen(f'{0}/makeblastdb -in {knownseqfile} -dbtype prot -hash_index -out {knownseqfile}'.format(
    #     config.C['ncbi_path'],))


def check(knownseqfile):
    #Check if an NCBI db of the knownseqfile already exists
    query = knownseqfile

    if os.path.isfile(f'{knownseqfile}.phd') :
        print(f" [BLAST] NCBI database for {knownseqfile} exists already !\n")
    else:
        print(f" [BLAST] Creating an NCBI database for {knownseqfile} in {config.C['neuropeptide_family']}.\n")
        make(knownseqfile)
    return query


def blastp(blastdb, queryfile, outfile):
    print(' [BLAST] Running BLASTp...')
    blastp = os.popen(f'blastp -db {blastdb} -query {queryfile} -outfmt 10 -out {outfile}')
    # blastp = os.popen('{0}/blastp -db {1} -query {2} -outfmt 10 -out {3}'.format(
    #     config.C['ncbi_path'],blastdb,queryfile,outfile))
    print('  [BLAST] BLASTp complete.')

    return blastp


def file(file):
    First_line = True
    F = []
    with open(file) as f:
        for l in f:
            ll = l.strip().split(',')
            if First_line is True:
                First_line = False
                continue
            F.append(ll)

    return F


def parse(blastoutput):
    #Parsing the blast outfmt 10 file"""
    B = []
    blastout = file(blastoutput)
    for b in blastout:
        #[hit,known sequence,PID,length,evalue]
        B.append([int(b[0]),b[1],b[2],b[3],b[10]])

    #Create list of unique query id's"""
    query = set([i[0] for i in B])
    U = []
    for q in query:
        #Loop through blastp output and group query ID's"""
        b = [i for i in B if i[0] == q]
        
        #For each query ID, sort group according to evalue"""
        b.sort(key=lambda i:float(i[4]))

        #For each group, keep only the hit with the lowest evalue"""
        U.append(b[0])
    U.sort()
    
    return U
