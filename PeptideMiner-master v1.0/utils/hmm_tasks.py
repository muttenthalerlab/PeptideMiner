import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------
def run_hmmsearch(Output,HmmFile,Query):
# -----------------------------------------------------------------------
    cmd = f"hmmsearch --tblout {Output}.tbl {HmmFile} {Query}" 
    
    logger.info(f" [HMM Search] -> {Output}.tbl")
    p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
    ret = p.stdout
    
    # run = sub.Popen(_cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
    # stdout,stderr = run.communicate()
    #return [stdout,stderr]
    
    return [ret]

# -----------------------------------------------------------------------
def addsequence_hmmsearch(hmmsearch_dir,hmmsearch_file,dict_Fasta):
# -----------------------------------------------------------------------
    logger.info(f" [HMM Search] -> {hmmsearch_file}")
    
    csv_header =["ID", "accession", "query_name", "accession", "full_E-value", "full_score", "full_bias", 
                "dom_E-value", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", 
                "inc", "desc_target", "sequence"]
    
    with open(f"{os.path.join(hmmsearch_dir,hmmsearch_file)}.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_header)
        
        with open(f"{os.path.join(hmmsearch_dir,hmmsearch_file)}.tbl", "r") as tblfile:
            for line in tblfile.readlines():
                if not line.startswith("#"):
                    ll = line.replace(',','').strip().split()
                    for ID in dict_Fasta:
                        if str(ID).startswith(str(ll[0])):
                            ll.append(dict_Fasta[ID])
                            writer.writerow(ll)

    return(f"{hmmsearch_file}.csv")


# -----------------------------------------------------------------------
def filter_hmmsearch(hmmsearch_dir,hmmsearch_csv,hmm_name, transcriptome_name):
# -----------------------------------------------------------------------
    with open(os.path.join(hmmsearch_dir,hmmsearch_csv), "r") as f:    
        reader = csv.DictReader(f)
        hmm_search = list(reader)

    hmm_search_sort = sorted(hmm_search, key=lambda d: d['full_E-value'])
    unique_id = set([i['ID'] for i in hmm_search_sort])
    
    lowest_hmm_search = {}
    for hmm in hmm_search_sort:
        if hmm['ID'] not in lowest_hmm_search:
            lowest_hmm_search[hmm['ID']] = [hmm['ID'],hmm_name, transcriptome_name, hmm['full_E-value'],hmm['sequence']]
    return(lowest_hmm_search)
        
