import os, csv
import datetime

import logging
logger = logging.getLogger(__name__)

from lib.database import sql_connector
from lib.db_tasks import upload_known_peptides, upload_cds, get_seqreads, get_summary_familyname

# ====================================================================================================
class PeptideMiner():
# ====================================================================================================

    def __init__(self, prgArgs):

        # Run Parameters
        self.cds_min_length = int(prgArgs.cds_min_length)
        self.signalp_cutoff = float(prgArgs.signalp_cutoff)
        self.signalp_min_length = int(prgArgs.signalp_min_length)

        self.mature_evalue_cutoff = float(prgArgs.mature_evalue_cutoff)
        self.mature_min_length = int(prgArgs.mature_min_length)
        self.mature_max_length = int(prgArgs.mature_max_length)

        # Data Folders
        self.data_dir = prgArgs.datadir
        self.known_pep_dir = os.path.join(self.data_dir,'01-known_seq')
        self.hmm_dir = os.path.join(self.data_dir,'02-pHMM')
        self.database_file = 'sqlite.db'
        self.database_sql = 'PeptideMiner.sql'

        # Work Folders
        self.work_dir = prgArgs.workdir
        self.hmmsearch_dir = os.path.join(self.work_dir,'01-hmmsearch')
        self.pipeline_dir = os.path.join(self.work_dir,'02-pipeline')
        
        # Programs
        self.hmmsearch_path = 'hmmsearch'
        self.fasta36_path = 'fasta36'
        self.signalp_path = prgArgs.signalp_path

        # Pipeline
        self.family_name = prgArgs.peptide_family
        self.query_dir = prgArgs.querydir

        self.pipeline_filename = {
            '01': {'filename': '01_hmmsearch',    },
            '02': {'filename': '02_hmmsearch_seq',    },
            '03': {'filename': '03_cds',          },
            '04': {'filename': '04_signalp_seq',  },
            '05': {'filename': '05_mature_seq', },
            '06': {'filename': '06_mature_rmduplicate_seq', },
            '07': {'filename': '07_blastp', },
            '08': {'filename': '08_summary', },
        }
 
        # Pipeline Properties
        self.hmm_lst = []
        self.hmm_seq_lst = []
        self.knownpep_lst = []
        self.hmm_search_files = []
        self.hmm_query_files = []
        self.cds_lst = []
        self.maturepep_lst = []
        self.matureseq_lst = []
        self.blastp_annotations = []
        #self.seq_id_dict = {}

        # Initialise Working Folders
        if not os.path.exists(self.hmmsearch_dir):
            os.makedirs(self.hmmsearch_dir)
        if not os.path.exists(self.pipeline_dir):
            os.makedirs(self.pipeline_dir)

        # SQL Database
        self.sql_database_file = os.path.join(self.data_dir,self.database_file)
        self.sql_create = os.path.join(self.data_dir,self.database_sql)
        self.db = sql_connector(data_file= self.sql_database_file, sql_create= self.sql_create)

    # ----------------------------------------------------------
    @staticmethod
    def read_fasta_file(FastA_File):
        """
         Reads FastA file into Dict with >line as key and sequence as value
        """
        dict_Seq = {}
        kSeq = None
        if os.path.isfile(FastA_File):
            with open(FastA_File,encoding='latin-1') as f:
                for line in f:
                    if line.startswith('>'):
                        kSeq = line[1:]
                        dict_Seq[kSeq] = "" 
                    else:
                        dict_Seq[kSeq] += line.strip()
        return(dict_Seq)

    # ----------------------------------------------------------
    @staticmethod
    def write_fasta_file(Fasta_File,Fasta_Dict):
        """
         Writes FastA file from Dict with >line as key and sequence as value
        """
        with open(Fasta_File,'w') as f:
            for key in Fasta_Dict:
                #print(f" [write_fasta_file] - [{key}]")
                f.write(f">{key}\n{Fasta_Dict[key]}\n")

    # ----------------------------------------------------------
    def check_paths(self):
        Error_Dict = {}
        if not os.path.isfile(self.signalp_path):
            Error_Dict['SignalP'] = f"Not Found {self.signalp_path}"
        return(Error_Dict)

# ====================================================================================================
def read_known_peptides(PM):
# ====================================================================================================
# Step 1.1 - Read known peptide sequences
#
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
#
    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['03']['filename']}_sequences.csv"
    csv_header=['cds_id','seq_id','n_cds','cds']

    # All Sequences or just from the specific hmm_id
    seq_id_dict = get_seqreads(PM.db) 

    _cds_lst = []
    for seq_id in seq_id_dict:
        n_cds = 0
        for seq_seg in seq_id['precursor'].split('*'):
            seq_M = seq_seg
            
            if len(seq_seg) >= PM.cds_min_length and 'M' in seq_seg:
                seq_M = seq_seg[seq_seg.index('M'):]
                if len(seq_M) >= PM.cds_min_length:
                    n_cds += 1                
                    _cds_lst.append({'seq_id':seq_id['id'],'n_cds':n_cds,'cds':seq_M})
        if n_cds == 0:
            # if no cds use original precursor
           _cds_lst.append({'seq_id':seq_id['id'],'n_cds':0,'cds':seq_id['precursor']})
    
    # Upload and make unique by cds_id
    PM.cds_lst = []        
    _cds_ids = []
    for cds in _cds_lst:
        cds_id = upload_cds(PM.db,cds['cds'],cds['seq_id'])
        if cds_id not in _cds_ids:
            cds['cds_id'] = cds_id
            PM.cds_lst.append(cds)
            _cds_ids.append(cds_id)
                   
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

    SumTime= datetime.datetime.now()

    csv_dir = PM.pipeline_dir
    csv_filename = f"{PM.pipeline_filename['08']['filename']}_{PM.family_name}_sequences.csv"
    txt_filename = f"{PM.pipeline_filename['08']['filename']}_{PM.family_name}.txt"

    sum_data = get_summary_familyname(PM.db,PM.family_name)
    print(f" [Summary] {len(sum_data)}")
    if sum_data:
        csv_header = list(sum_data[0].keys())
        # Write CSV file

        logger.info(f" [BlastP] Annotations -> {csv_filename} ({len(PM.blastp_annotations)})")
        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
            csvwriter.writeheader()
            for s in sum_data:
                csvwriter.writerow(s)

        _set_peptideminer_hits = set(l['hit id'] for l in sum_data)
        _set_querydb = set(l['hit query DB'] for l in sum_data)
        _set_phmm = set(l['pHMM name'] for l in sum_data)
        _set_unique_matseq = set(l['hit mature sequence'] for l in sum_data)
        _set_unique_pre = set(l['hit CDS'] for l in sum_data)
        _set_hmm = set(l['hmm_name'] for l in PM.hmm_lst)

        with open(os.path.join(csv_dir,txt_filename),'w') as out:
            out.write(f"Summary PeptideMiner peptide search\n")
            out.write(f"Date: {SumTime.strftime('%d/%m/%Y')}\n")
            out.write("\n")
            out.write(f"Number of profile-HMMs used: {len(_set_hmm)}\n")
            out.write(f"\t{','.join(_set_hmm)}\n")
            out.write(f"Number of databases searched: {len(PM.hmm_query_files)}\n")
            out.write(f"\t{','.join(PM.hmm_query_files)}\n")
            out.write("\n")
            out.write(f"Output:\n")
            out.write(f"hmmsearch hits: {len(PM.hmm_seq_lst)}\n") #
            out.write(f"PeptideMiner hits: {len(_set_peptideminer_hits)}\n")
            out.write(f"\tNumber of unique CDS: {len(_set_unique_pre)}\n")
            out.write(f"\tNumber of unique mature peptides: {len(_set_unique_matseq)}\n")
    else:
        logger.error(f" [BlastP] NO Hits found !!")