#
import subprocess

import logging
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------
def upload_known_peptides(DB,family_name,species,accession,seq_name,seq, verbose=0):
# -----------------------------------------------------------------------

    # Check id for family name
    family_id = DB.onevalue('neuropeptide_family','id',{'name':family_name})
    if family_id is None:
        DB.insert('neuropeptide_family',{'name':family_name})
        family_id = DB.onevalue('id','neuropeptide_family',{'name':family_name})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [neuropeptide_family] new {family_name} ({family_id})")

    # Add new Sequence if not exists
    peptide_id = DB.onevalue('known_NP','id',{'familyid':family_id,'accession':accession,'name':seq_name})
    if peptide_id is None:
        DB.insert('known_NP',{'name':seq_name,'familyid':family_id,'species':species,'sequence':seq,'accession':accession})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [known_NP] new {seq_name} for {species}")

    return(family_id,peptide_id)


# -----------------------------------------------------------------------
def upload_hmmsearch(DB,hmm_name,transcriptome_name,id,evalue,sequence, verbose=0):
# -----------------------------------------------------------------------

    hmm_id = DB.onevalue('hmm','id',{'name':hmm_name})

    if hmm_id is None:
        DB.insert('hmm',{'name':hmm_name})
        hmm_id = DB.onevalue('id','hmm',{'name':hmm_name})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [hmm] new {hmm_name} ({hmm_id})")

    seqreads_id = DB.onevalue('seqreads','id',{'name':id,'hmmid':hmm_id})
    if seqreads_id is None:
        DB.insert('seqreads',{'name':id,'hmmid':hmm_id,'transcriptome':transcriptome_name,'evalue':evalue,'signalseq_length':0,'precursor':sequence})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [seqreads] new {id} for {hmm_name} {transcriptome_name}")
        # readname = DB.onevalue('name','seqreads',{'name':id})
        # hmmid_seqread = DB.onevalue('hmmid','seqreads',{'hmmid':hmm_id})

    return(hmm_id,seqreads_id)


# -----------------------------------------------------------------------
def upload_cds(DB,cds,seq_id,verbose=0):
# -----------------------------------------------------------------------

    cds_id = DB.onevalue('cds','id',{'sequence':cds})
        
    if cds_id is None:
        DB.insert('cds',{'sequence':cds})
        cds_id = DB.onevalue('cds','id',{'sequence':cds})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [cds] new {cds_id} {len(cds)} AA ")

    DB.update('seqreads',{'cds_id':cds_id},f'id={seq_id}')
    return cds_id



# -----------------------------------------------------------------------
def get_hmmid(DB,hmm_name):
# -----------------------------------------------------------------------
    for res in DB.selectall(('id'),'hmm',(f"name='{hmm_name}'")):
        return str(res[0])
    
# -----------------------------------------------------------------------
def get_seqreads(DB):
# -----------------------------------------------------------------------    
    lst_Dict = []
    for res in DB.selectall('seqreads s',('id','precursor','hmmid')):
        lst_Dict.append({'id':res[0],'precursor':res[1],'hmmid': res[2]})
    return lst_Dict
