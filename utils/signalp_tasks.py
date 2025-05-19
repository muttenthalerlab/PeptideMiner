import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

SIGNALP_HEADER = ['name','Cmax','CMax_pos','Ymax','Ymax_pos','Smax','Smax_pos','Smean','Dscore','SP','Dmaxcut','Networks-used']

# -----------------------------------------------------------------------
def run_signalp(seq_id,sequence,cutoff,signal_path,temp_path='/tmp'):
# -----------------------------------------------------------------------
    # Create sequence fasta file
    seq_file = os.path.join(temp_path,f"{seq_id}.fasta")
    with open(seq_file,'w') as tmp:
        tmp.write(f">{seq_id}\n{sequence}")
    
    #out_file = os.path.join(temp_path,f"{seq_id}.out")

    # SignalP 4.1  : eukariotic 
    cmd = f"{signal_path} -t euk -M -U {cutoff} -u {cutoff} {seq_file} "
    p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
    ret = p.stdout
    os.remove(seq_file)

    signalp_dict = dict(zip(SIGNALP_HEADER,ret.splitlines()[2].split()))

    logger.info(f" [SignalP] -> {seq_id} {signalp_dict['SP']} at {signalp_dict['CMax_pos']}")
    return(signalp_dict)




