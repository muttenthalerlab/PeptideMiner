import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

from lib.db_tasks import upload_annotations

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
    qry_outlst = []
    if os.path.isfile(QueryOut):
        with open(QueryOut,'r') as f:
            csvreader = csv.DictReader(f, fieldnames= BLASTP_QRY_HEADER)
            for row in csvreader:
                qry_outlst.append(row)
    return(qry_outlst)


# ====================================================================================================
def run_blast(PM, Overwrite=False):
# ====================================================================================================

    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['07']['filename']}_annotation.csv"

    blast_filename = f"{PM.pipeline_filename['07']['filename']}_db"

    qry_fasta = f"{PM.pipeline_filename['06']['filename']}.fna"  
    qry_out = f"{PM.pipeline_filename['07']['filename']}_query.csv"
    
    # Write Fasta file
    _fasta = {}
    for s in PM.knownpep_lst:
        _name = ":".join([str(s['family_id']),str(s['peptide_id']),s['species'],s['name'],s['accession']])
        _fasta[_name.strip()] = s['sequence']
    logger.info(f" [BlastP] MatureSeq: -> {blast_filename}.fna ({len(_fasta)})")
    PM.write_fasta_file(f"{os.path.join(csv_dir,blast_filename)}.fna",_fasta)

    # Make BlastP database
    make_blastp_db(os.path.join(csv_dir,blast_filename))

    # Run and Parse BlastP query for unique qry_names with lowest evalue
    run_blastp_query(os.path.join(csv_dir,blast_filename),
                     os.path.join(csv_dir,qry_fasta),
                     os.path.join(csv_dir,qry_out))
    blastp_out = parse_blastp_query(os.path.join(csv_dir,qry_out))
    logger.info(f" [BlastP] BastP Out {len(blastp_out)} ")
    unq_qry_names = set([b['qry_name'] for b in blastp_out])
    
    logger.info(f" [BlastP] UnqQryNames {len(unq_qry_names)} ")

    PM.blastp_annotations = []
    for qn in unq_qry_names:
        _q = [i for i in blastp_out if i['qry_name'] == qn]
        _q.sort(key=lambda i:float(i['evalue']))
        PM.blastp_annotations.append(_q[0])
    
    # Upload Annotation
    for a in PM.blastp_annotations:
        _nid = a['qry_name'].split(':')[0]
        _kid = a['subject_name'].split(':')[0]
        upload_annotations(PM.db,_nid,_kid,float(a['pct_identity']),float(a['evalue']),int(a['length']))
    logger.info(f" [BlastP] Annotations {len(PM.blastp_annotations)} uploded")

    # Write CSV file
    logger.info(f" [BlastP] Annotations -> {csv_filename} ({len(PM.blastp_annotations)})")
    with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
        csvwriter = csv.DictWriter(f, fieldnames=BLASTP_QRY_HEADER)                
        csvwriter.writeheader()
        for a in PM.blastp_annotations:
            csvwriter.writerow(a)
