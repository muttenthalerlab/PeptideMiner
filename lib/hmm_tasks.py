import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

from lib.db_tasks import upload_hmmsearch

HMM_SEARCH_HEADER =["ID", "accession", "query_name", "accession", "full_E-value", "full_score", "full_bias", 
                    "dom_E-value", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", 
                    "inc", "desc_target", "sequence"]


# -----------------------------------------------------------------------
def run_hmmsearch(Output,HmmFile,Query,HMMSearch_Path='hmmsearch'):
# -----------------------------------------------------------------------
    cmd = f"{HMMSearch_Path} --tblout {Output}.tbl {HmmFile} {Query}" 
    
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
        

# ====================================================================================================
def hmmsearch(PM, Overwrite=False):
# ====================================================================================================
# Step 1.2 HMM search
#
    # Output Summary        
    hmm_search_outfile = os.path.join(PM.hmmsearch_dir,f"{PM.pipeline_filename['01']['filename']}.txt") 
    hmm_search_out = []
    
    # Check if hmmsearch has been run before
    PM.hmm_search_files = [f for f in os.listdir(PM.hmmsearch_dir) if f.endswith('csv')]
    
    if len(PM.hmm_search_files)==0 or Overwrite:
        PM.hmm_search_files = []
        
        #Get a list of the profile-HMMs available for the desired family
        HMM_Files = [h for h in os.listdir(PM.hmm_dir) if h.startswith(PM.family_name)]
        logger.info(f" [HMM Search] {PM.family_name} - HMM Files : {HMM_Files}")
        
        #For each Query (fna) file run hmmsearch
        PM.hmm_query_files = [q for q in os.listdir(PM.query_dir) if q.endswith('fna')]
        for qry_file in PM.hmm_query_files:
            
            dict_Fasta = PM.read_fasta_file(os.path.join(PM.query_dir,qry_file))
            logger.info(f" [HMM Search] {PM.family_name} - HMM Query : {qry_file} with {len(dict_Fasta)} sequences")

            for hmm in HMM_Files:
                hmm_out = f"{hmm.split('.')[0]}.{qry_file.replace('.fna','')}"
                
                # Run HMM Search 
                _ret = run_hmmsearch(os.path.join(PM.hmmsearch_dir,hmm_out),
                                     os.path.join(PM.hmm_dir,hmm),
                                     os.path.join(PM.query_dir,qry_file),
                                     HMMSearch_Path=PM.hmmsearch_path)
                hmm_search_out.append(_ret)
                                            
                # Add Seqeunce infor to HMM Search Output CSV file
                addsequence_hmmsearch(PM.hmmsearch_dir,hmm_out,dict_Fasta)
                PM.hmm_search_files.append(f"{hmm_out}.csv")
        
        # Output logfiles
        with open(hmm_search_outfile,'w') as out:
            for line in hmm_search_out:
                out.writelines(line)
    else:
        logger.info(f" [HMM Search] Found {len(PM.hmm_search_files)} HMM Searches in {PM.hmmsearch_dir}")

# ====================================================================================================
def read_hmmsearch(PM, Overwrite=False):
# ====================================================================================================
# Step 1.2 Read HMMsearch
#
    
    if len(PM.hmm_search_files)>0:
        logger.info(f" [HMM Search] Reading {len(PM.hmm_search_files)} HMM Searches ")
        
        hmm_id_dict = {}
        seq_id_dict = {}
                    
        for hmm in  PM.hmm_search_files:
            _hmm_lst = hmm.split('.')                
            hmm_name = _hmm_lst[0]
            transcriptome_name =_hmm_lst[1]
            
            hmm_lstdict = filter_hmmsearch(PM.hmmsearch_dir,hmm,hmm_name,transcriptome_name)
            for id in hmm_lstdict:
                evalue = hmm_lstdict[id][3]
                sequence = hmm_lstdict[id][4]
                hmm_id, seq_id = upload_hmmsearch(PM.db,hmm_name,transcriptome_name,id,evalue,sequence)

                hmm_id_dict[hmm_id] = [hmm_id,hmm_name,transcriptome_name]
                seq_id_dict[seq_id] = [seq_id,sequence,hmm_id,hmm_name]

    # Convert dict_lst into lst_dict
    PM.hmm_lst = []
    for k in hmm_id_dict:
        PM.hmm_lst.append({'hmm_id':hmm_id_dict[k][0],'hmm_name':hmm_id_dict[k][1],'transcriptome':hmm_id_dict[k][2]})

    PM.hmm_seq_lst = []
    for k in seq_id_dict:
        PM.hmm_seq_lst.append({'seq_id':seq_id_dict[k][0],'sequence':seq_id_dict[k][1],'hmm_id':seq_id_dict[k][2],'hmm_name':seq_id_dict[k][3]})

    logger.info(f" [HMM Search] HMM's: {len(hmm_id_dict)} uploaded")
    logger.info(f" [HMM Search] Seq's: {len(seq_id_dict)} uploaded")
        
    # Pipeline CSV Log
    csv_filename = f"{PM.pipeline_filename['01']['filename']}.csv" 
    with open(os.path.join(PM.pipeline_dir,csv_filename),'w',newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csv_header=['hmm_id','hmm_name','transcriptome']
        csvwriter.writerow(csv_header)
        for key in hmm_id_dict:
            csvwriter.writerow(hmm_id_dict[key])
    logger.info(f" [HMM Search] HMM's -> {csv_filename} ({len(hmm_id_dict)})")

    csv_filename = f"{PM.pipeline_filename['01']['filename']}_sequences.csv"
    with open(os.path.join(PM.pipeline_dir,csv_filename),'w',newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csv_header=['seq_id','sequence','hmm_id','hmm_name']
        csvwriter.writerow(csv_header)
        for key in seq_id_dict:
            csvwriter.writerow(seq_id_dict[key])
    logger.info(f" [HMM Search] Seq's -> {csv_filename} ({len(seq_id_dict)})")
