import os
from libfunc import mysqlpop, config


"""
Main script support module.
Part of the SQLite DB set up.
Populates known_seq and neuropeptide_family table.
BLAST noduplicates against known_seq table.
"""

def known(file):
    S = []
    species = ''
    accession = ''
    name = ''
    seq = ''
    with open(file, "r", encode='latin-1') as f:
        for l in f:
            if l.startswith('>'):
                S.append([species,accession,name,seq])
                species = l.strip().split('[')[-1].replace(']','').split('(')[0]
                accession = l.split('_')[0][1:]
                name = ' '.join(l.split('[')[0].split('_')[1:])
            else:
                seq = l.strip()
        S.append([species,accession,name,seq])
        S.pop(0)
    return S

def known_pop():
    count = 0
    knownseqfile = ''
    knownseqpath = f"{config.C['path']}/data/01-known_seq/{config.C['neuropeptide_family']}"

    for k in os.listdir(knownseqpath):
        if k.endswith('.fna'):
            knownseqfile = f"{knownseqpath}/{k}"

    if knownseqfile == '':
        print(f"No known peptide file in directory {knownseqpath}")
        print(f"The known peptide file must be a fasta file and have the extension .fna.")
        exit()
    else:
        print(f"Reading known peptides file {knownseqfile}")
        for s in known(knownseqfile):
            c = mysqlpop.knownseq(s,config.C['neuropeptide_family'])
            count += c
    print("The Known_NP table in the MySQL database is populated")
    print(f"{count} known peptides were added to the SQLite Known_NP table.\n") 
