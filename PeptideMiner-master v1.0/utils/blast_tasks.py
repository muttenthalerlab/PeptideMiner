import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

BLASTP_QRY_HEADER = ['qry_name','subject_name',
                     'pct_identity','length','n_mismatch','n_gap',
                     'q_start','q_end','s_start','s_end',
                     'evalue','hitscore']

# -----------------------------------------------------------------------
def make_blastp_db(FileName, Overwrite=False, verbose=0):
# -----------------------------------------------------------------------
    ret = None
    if not os.path.isfile(f'{FileName}.phd') or Overwrite:
        # BlastP  
        cmd = f"makeblastdb -in {FileName}.fna -dbtype prot -hash_index -out {FileName} "
        logger.info(f" [BlastP] DB(Prot) -> {FileName} ")
        p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
        ret = p.stdout
    return(ret)

# -----------------------------------------------------------------------
def run_blastp_query(BlastDB, QueryFasta, QueryOut):
# -----------------------------------------------------------------------
    ret = None
    cmd = f"blastp -db {BlastDB} -query {QueryFasta} -outfmt 10 -out {QueryOut}"
    logger.info(f" [BlastP] Search -> {QueryOut} ")
    p = subprocess.run(cmd,shell=True,capture_output=False, text=True)
    ret = p.stdout
    return(ret)

# -----------------------------------------------------------------------
def parse_blastp_query(QueryOut):
# -----------------------------------------------------------------------
    # BlastP_Header = ['qry_name','subject_name',
    #                  'pct_identity','length','n_mismatch','n_gap',
    #                  'q_start','q_end','s_start','s_end',
    #                  'evalue','hitscore']
    
    qry_outlst = []
    if os.path.isfile(QueryOut):
        with open(QueryOut,'r') as f:
            csvreader = csv.DictReader(f, fieldnames= BLASTP_QRY_HEADER)
            for row in csvreader:
                qry_outlst.append(row)
    return(qry_outlst)
