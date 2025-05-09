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
        family_id = DB.onevalue('neuropeptide_family','id',{'name':family_name})
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
    # Check id for HMM name
    hmm_id = DB.onevalue('hmm','id',{'name':hmm_name})
    if hmm_id is None:
        DB.insert('hmm',{'name':hmm_name})
        hmm_id = DB.onevalue('hmm','id',{'name':hmm_name})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [hmm] new {hmm_name} ({hmm_id})")

    # Check id for Sequences
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
def upload_matureseq(DB,cds_id,seq, verbose=0):
# -----------------------------------------------------------------------        
    mature_id = DB.onevalue('mature','id',{'cds_id':cds_id,'matseq':seq})
    if mature_id is None:
        DB.insert('mature',{'cds_id':cds_id,'matseq':seq})
        mature_id = DB.onevalue('mature','id',{'cds_id':cds_id,'matseq':seq})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [mature] new {mature_id} {len(seq)} AA ")
    return mature_id

# -----------------------------------------------------------------------
def upload_noduplicates(DB,cds_id,seq, verbose=0):
# -----------------------------------------------------------------------        
    transcriptome_cdsid = DB.onevalue('seqreads','transcriptome',{'cds_id':cds_id})
    hmmid_cdsid = DB.onevalue('seqreads','hmmid',{'cds_id':cds_id})

    nodup_id = DB.onevalue('noduplicates','id',{'transcriptome':transcriptome_cdsid,'matseq':seq})
    if nodup_id is None:
        DB.insert('noduplicates',{'hmm_id':hmmid_cdsid,'transcriptome':transcriptome_cdsid,'matseq':seq})        
        nodup_id = DB.onevalue('noduplicates','id',{'matseq':seq})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [noduplicates] new {nodup_id} {len(seq)} AA ")

    # Update 'mature' table with new id from 'noduplicates' table
    DB.update('mature',{'noduplicates_id':nodup_id},f"matseq='{seq}'")
    return(nodup_id)

# -----------------------------------------------------------------------
def upload_annotations(DB,novel_id,known_id,pct_identity,evalue,length, verbose=0):
# -----------------------------------------------------------------------
    annotated_id = DB.onevalue('annotated','id',{'novel_id':novel_id,'knownNP_id':known_id})
    if annotated_id is None:
        DB.insert('annotated',{'novel_id':novel_id,'knownNP_id':known_id,'pid':pct_identity,'evalue':evalue,'length_alignment':length})
        annotated_id = DB.onevalue('annotated','id',{'novel_id':novel_id,'knownNP_id':known_id})
        if verbose > 0:
            logger.info(f" [SQLlite] Table [annotated] new {annotated_id} {length} AA ")
    return(annotated_id)

# -----------------------------------------------------------------------
def update_seqreads_signalp(DB,cds_id,pos):
# -----------------------------------------------------------------------
    DB.update('seqreads',{'signalseq_length':pos},f'cds_id={cds_id}')

# -----------------------------------------------------------------------
def get_hmmid(DB,hmm_name):
# -----------------------------------------------------------------------
    for res in DB.selectall('hmm',('id'),(f"name='{hmm_name}'")):
        return str(res[0])

# -----------------------------------------------------------------------
def get_noduplicates(DB,hmm_id):
# -----------------------------------------------------------------------    
    out = []
    for res in DB.selectall('noduplicates',('id','hmm_id','transcriptome','matseq'),(f"hmm_id={hmm_id}")):
        out.append({'id':res[0], 'hmm_id': res[1], 'transcriptome': res[2], 'matseq': res[3]})
    return out
    
# -----------------------------------------------------------------------
def get_seqreads(DB):
# -----------------------------------------------------------------------    
    lst_Dict = []
    for res in DB.selectall('seqreads s',('id','precursor','hmmid')):
        lst_Dict.append({'id':res[0],'precursor':res[1],'hmmid': res[2]})
    return lst_Dict
