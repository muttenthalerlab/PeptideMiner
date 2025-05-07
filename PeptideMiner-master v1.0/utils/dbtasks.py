#
import subprocess

import logging
logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------
def upload_known_peptides(DB,family_name,species,accession,seq_name,seq):
# -----------------------------------------------------------------------

    # Check id for family name
    family_id = DB.onevalue('id','neuropeptide_family',{'name':family_name})
    if family_id is None:
      DB.insert('neuropeptide_family',{'name':family_name})
      family_id = DB.onevalue('id','neuropeptide_family',{'name':family_name})
      logger.error(f" [SQLlite] Table [neuropeptide_family] new {family_name} ({family_id})")

    # Add new Sequence if not exists
    peptide_id = DB.onevalue('id','known_NP',{'familyid':family_id,'accession':accession,'name':seq_name})
    if peptide_id is None:
      DB.insert('known_NP',{'name':seq_name,'familyid':family_id,'species':species,'sequence':seq,'accession':accession})
      logger.error(f" [SQLlite] Table [known_NP] new {seq_name} for {species}")
      return(1)
    return(0)


# -----------------------------------------------------------------------
def run_hmmsearch(OutPut,HmmFile,Query):
# -----------------------------------------------------------------------
    cmd = f"hmmsearch --tblout {OutPut}.tbl {HmmFile} {Query}" 
    
    logger.info(f"[HMM Search] {Query} {HmmFile} -> {OutPut}")
    p = subprocess.run(cmd,shell=True,capture_output=True, text=True)
    ret = p.stdout
    
    # run = sub.Popen(_cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
    # stdout,stderr = run.communicate()
    #return [stdout,stderr]
    
    return [ret]