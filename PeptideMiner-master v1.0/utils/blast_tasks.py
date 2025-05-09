import os
import subprocess

import logging
logger = logging.getLogger(__name__)

def make_blastp_db(FileName, Overwrite=False):
    ret = None
    if not os.path.isfile(f'{FileName}.phd') or Overwrite:
        # BlastP  
        cmd = f"makeblastdb -in {FileName} -dbtype prot -hash_index -out {FileName} "
        logger.info(f" [BlastP] DB(Prot) -> {FileName} ")
        p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
        ret = p.stdout
    return(ret)