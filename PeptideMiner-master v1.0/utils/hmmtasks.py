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

    # csv_file = open(f"{os.path.join(hmmsearch_dir,hmmsearch_file)}.csv", "w")
    # csv_file.writelines(','.join(csv_header))         
    # for line in open(f"{os.path.join(hmmsearch_dir,hmmsearch_file)}.tbl", "r").readlines():
    #     if not line.startswith("#"):
    #         ll = line.replace(',','').strip().split()
    #         for ID in dict_Fasta:
    #             #If ID in .tbl matches an ID in fastaDict
    #             if str(ID).startswith(str(ll[0])):
    #                 csv_line = ",".join(ll[:18]) + "," + " ".join(ll[18:]) + "," + str(dict_Fasta[ID])
    #                 csv_file.writelines(csv_line)
    # csv_file.close()
    return(f"{hmmsearch_file}.csv")


# -----------------------------------------------------------------------
def filter_hmmsearch(hmmsearch_dir,hmmsearch_csv,hmm_name, transcriptome_name):
# -----------------------------------------------------------------------
    with open(os.path.join(hmmsearch_dir,hmmsearch_csv), "r") as f:    
        reader = csv.DictReader(f)
        hmm_search = list(reader)

    print(hmm_search)    
    # E=[]
    # csv_file = str(input_file).split('/')[-1]
    # hmm = str(csv_file).split('.')[0]
    # transcriptome_name = '.'.join(str(csv_file).split('.')[1:-1])

    # first_line = True
    # with open(input_file) as f:
    #     for l in f:
    #         fields = l.strip().split(',')
    #         if first_line:
    #             first_line = False 
    #             continue
    #         E.append([fields[0],hmm,transcriptome_name,fields[4],fields[19]])

    # #Isolate readnames in the table and keep only the unique ones
    # nodup = set([i[0] for i in E])
    
    # S = []
    # for r in nodup:
    #     #Isolate lines where the r == readname
    #     b = [i for i in E if i[0] == r]
        
    #     #Sort list accroding to evalue
    #     b.sort(key=lambda i:float(i[3]))
    #     S.append(b[0])
    # return S
