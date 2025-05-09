import os
import subprocess

import logging
logger = logging.getLogger(__name__)

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
def run_blastp_search(BlastDB, QueryFasta, QueryOut):
# -----------------------------------------------------------------------
    ret = None
    cmd = f"blastp -db {BlastDB} -query {QueryFasta} -outfmt 10 -out {QueryOut}"
    logger.info(f" [BlastP] Search -> {QueryOut} ")
    p = subprocess.run(cmd,shell=True,capture_output=False, text=True)
    ret = p.stdout
    return(ret)
