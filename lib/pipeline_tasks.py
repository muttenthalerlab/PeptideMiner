import os, csv

import logging
logger = logging.getLogger(__name__)

from lib.db_tasks import upload_known_peptides, upload_cds, get_seqreads, get_summary_familyname

# ====================================================================================================
def read_known_peptides(PM):
# ====================================================================================================
    PM.known_peptide = []
    
    pep_dir = os.path.join(PM.known_pep_dir,PM.family_name)
    for k in os.listdir(pep_dir):
        if k.endswith('.fna'):
            PM.known_peptide.append(os.path.join(pep_dir,k))

    for fna_file in PM.known_peptide:
        if os.path.isfile(fna_file):
            logger.info(f" [Known Peptides] {PM.family_name} - Reading {fna_file}")
            dict_Seq = PM.read_fasta_file(fna_file)
            PM.n_known_peptide = 0
            for k in dict_Seq:
                # Parsing Fast Header
                #>P13204_ANFB [Bos taurus]
                #>acession_name [specie]
                #
                species = k.strip().split('[')[-1].replace(']','').split('(')[0]
                accession = k.split('_')[0][1:]
                name = ' '.join(k.split('[')[0].split('_')[1:])

                family_id,peptide_id = upload_known_peptides(PM.db, PM.family_name,species,accession,name,dict_Seq[k])
                PM.n_known_peptide += 1
                PM.knownpep_lst.append({'family_id':family_id,'peptide_id':peptide_id,
                                            'species':species,'accession':accession,'name':name,
                                            'sequence':dict_Seq[k]})
    else:
        logger.error(f" [Known Peptides] {PM.family_name} - {PM.n_known_peptide} peptide uploaded")

# ====================================================================================================
def read_cds(PM, Overwrite=False):
# ====================================================================================================

    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['03']['filename']}_sequences.csv"
    csv_header=['cds_id','seq_id','n_cds','cds']

    # All Sequences or just from the specific hmm_id
    seq_id_dict = get_seqreads(PM.db) 

    PM.cds_lst = []
    for seq_id in seq_id_dict:
        
        for seq_seg in seq_id['precursor'].split('*'):
            n_cds = 0
            seq_M = seq_seg
            
            if len(seq_seg) >= PM.cds_min_length and 'M' in seq_seg:
                seq_M = seq_seg[seq_seg.index('M'):]
                if len(seq_M) >= PM.cds_min_length:
                    n_cds += 1                
                    PM.cds_lst.append({'seq_id':seq_id['id'],'n_cds':n_cds,'cds':seq_M})
            
    for cds in PM.cds_lst:
        cds_id = upload_cds(PM.db,cds['cds'],cds['seq_id'])
        cds['cds_id'] = cds_id   
    logger.info(f" [HMM Search] CDS: {len(PM.cds_lst)} uploaded")

    with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
        csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
        csvwriter.writeheader()
        for cds in PM.cds_lst:
            csvwriter.writerow(cds)
    logger.info(f" [HMM Search] CDS: -> {csv_filename} ({len(PM.cds_lst)} )")

# ====================================================================================================
def summary(PM, Overwrite=False):
# ====================================================================================================

    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['08']['filename']}_{PM.family_name}.csv"

    sum_data = get_summary_familyname(PM.db,PM.family_name)
    csv_header = list(sum_data[0].keys())
    # Write CSV file

    logger.info(f" [BlastP] Annotations -> {csv_filename} ({len(PM.blastp_annotations)})")
    with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
        csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
        csvwriter.writeheader()
        for s in sum_data:
            csvwriter.writerow(s)
