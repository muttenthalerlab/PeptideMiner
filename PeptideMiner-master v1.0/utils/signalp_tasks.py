import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------
def run_signalp(seq_id,sequence,cutoff,signal_path,temp_path='/tmp'):

    # Create sequence fasta file
    seq_file = os.path.join(temp_path,f"{seq_id}.fasta")
    with open(seq_file,'w') as tmp:
        tmp.write(f">{seq_id}\n{sequence}")
    
    logger.info(f" [SingalP] -> ")
    # SignalP 4.1  : eukariotic 
    cmd = f"{signal_path} -t euk -U {cutoff} -u {cutoff} {seq_file}"
    p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
    ret = p.stdout

    print(ret)
