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

    return(family_id,peptide_id)


# -----------------------------------------------------------------------
def upload_hmmsearch(DB,hmm_name,transcriptome_name,id,evalue,sequence):
# -----------------------------------------------------------------------

    hmm_id = DB.onevalue('id','hmm',{'name':hmm_name})

    if hmm_id is None:
        DB.insert('hmm',{'name':hmm_name})
        hmm_id = DB.onevalue('id','hmm',{'name':hmm_name})
        logger.error(f" [SQLlite] Table [hmm] new {hmm_name} ({hmm_id})")

    seqreads_id = DB.onevalue('id','seqreads',{'name':id,'hmmid':hmm_id})
    if seqreads_id is None:
        DB.insert('seqreads',{'name':id,'hmmid':hmm_id,'transcriptome':transcriptome_name,'evalue':evalue,'signalseq_length':0,'precursor':sequence})
        logger.error(f" [SQLlite] Table [seqreads] new {id} for {hmm_name} {transcriptome_name}")
        # readname = DB.onevalue('name','seqreads',{'name':id})
        # hmmid_seqread = DB.onevalue('hmmid','seqreads',{'hmmid':hmm_id})

    return(hmm_id,seqreads_id)


# -----------------------------------------------------------------------
def get_hmmid(DB,hmm_name):
# -----------------------------------------------------------------------
    for res in DB.selectall(('id'),'hmm',(f"name='{hmm_name}'")):
        return str(res[0])
    
# -----------------------------------------------------------------------
def get_seqreads(DB):
# -----------------------------------------------------------------------    
    lst_Dict = []
    for res in DB.selectall(('id','precursor','hmmid'),'seqreads s'):
        lst_Dict.append({'id':res[0],'precursor':res[1],'hmmid': res[2]})
    return lst_Dict
