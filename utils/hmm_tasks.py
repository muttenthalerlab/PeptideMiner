import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)


HMM_SEARCH_HEADER =["ID", "accession", "query_name", "accession", "full_E-value", "full_score", "full_bias", 
                    "dom_E-value", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", 
                    "inc", "desc_target", "sequence"]


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
    
    with open(f"{os.path.join(hmmsearch_dir,hmmsearch_file)}.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(HMM_SEARCH_HEADER)
        
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
        

# ---------------------------------------------------------
def hmmsearch(PM, Overwrite=False):

    # Output Summary        
    hmm_search_outfile = os.path.join(PM.hmmsearch_dir,f"{PM.pipeline_filename['01']}.txt") 
    hmm_search_out = []
    
    # Check if hmmsearch has been run before
    PM.hmm_search_files = [f for f in os.listdir(PM.hmmsearch_dir) if f.endswith('csv')]
    
    if len(PM.hmm_search_files)==0 or Overwrite:
        PM.hmm_search_files = []
        
        #Get a list of the profile-HMMs available for the desired family
        HMM_Files = [h for h in os.listdir(PM.hmm_dir) if h.startswith(PM.family_name)]
        logger.info(f" [HMM Search] {PM.family_name} - HMM Files : {HMM_Files}")
        
        #For each Query (fna) file run hmmsearch
        Query_Files = [q for q in os.listdir(PM.query) if q.endswith('fna')]
        for qry_file in Query_Files:
            
            dict_Fasta = PM.read_fasta_file(os.path.join(PM.query,qry_file))
            logger.info(f" [HMM Search] {PM.family_name} - HMM Query : {qry_file} with {len(dict_Fasta)} sequences")

            for hmm in HMM_Files:
                hmm_out = f"{hmm.split('.')[0]}.{qry_file.replace('.fna','')}"
                
                # Run HMM Search 
                _ret = run_hmmsearch(os.path.join(PM.hmmsearch_dir,hmm_out),os.path.join(PM.hmm_dir,hmm),os.path.join(PM.query,qry_file))
                hmm_search_out.append(_ret)
                                            
                # Add Seqeunce infor to HMM Search Output CSV file
                addsequence_hmmsearch(PM.hmmsearch_dir,hmm_out,dict_Fasta)
                PM.hmm_search_files.append(f"{hmm_out}.csv")
        
        # Output logfiles
        with open(hmm_search_outfile,'w') as out:
            for line in hmm_search_out:
                out.writelines(line)
    else:
        logger.info(f" [HMM Search] Exists {PM.hmm_search_files} in {PM.hmmsearch_dir}")
