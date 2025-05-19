import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

from utils.db_tasks import update_seqreads_signalp

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

# ====================================================================================================
def signalp_cds(PM,Overwrite=False):
# ====================================================================================================

    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['04']['filename']}.csv"
    csv_header=['cds_id','signalp_pos','mature_peptide']

    if os.path.isfile(os.path.join(csv_dir,csv_filename)):
        PM.maturepep_lst = []
        with open(os.path.join(csv_dir,csv_filename)) as f:
            csvreader = csv.DictReader(f)
            for row in csvreader:
                PM.maturepep_lst.append(row)
        logger.info(f" [SignalP] MatPeptides: {csv_filename} -> ({len(PM.maturepep_lst)} )")
    else:
        PM.maturepep_lst = []
        if len(PM.cds_lst) > 0:
            for cds in PM.cds_lst:
                signalp_pos = 0
                mature_seq = cds['cds']
                        
                signalp_dict = run_signalp(f"{cds['cds_id']}_{cds['n_cds']:02d}",cds['cds'],PM.signalp_cutoff,PM.signalp_path)

                if signalp_dict['SP'] == 'Y':
                    signalp_pos = int(signalp_dict['CMax_pos'])
                    mature_seq = cds['cds'][signalp_pos:]
                    if len(mature_seq) >= PM.signalp_min_length:
                        PM.maturepep_lst.append({'cds_id':cds['cds_id'],'signalp_pos':signalp_pos,'mature_peptide':mature_seq})
                    else:
                        # ???? Not sure
                        PM.maturepep_lst.append({'cds_id':cds['cds_id'],'signalp_pos':signalp_pos,'mature_peptide':cds['cds']})

                update_seqreads_signalp(PM.db,cds['cds_id'],signalp_pos)

        
        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
            csvwriter.writeheader()
            for mpep in PM.maturepep_lst:
                csvwriter.writerow(mpep)
        logger.info(f" [SignalP] MatPeptides: -> {csv_filename} ({len(PM.maturepep_lst)} )")


