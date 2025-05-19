import os, sys
import datetime
import configargparse
from pathlib import Path

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "PeptideMiner"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(message)s",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)

from utils.database import sql_connector
from utils.pipeline_tasks import (read_known_peptides)
from utils.hmm_tasks import hmmsearch

# --------------------------------------------------------------------------------------
class PeptideMiner():
# --------------------------------------------------------------------------------------

    def __init__(self, prgArgs):

        # Run Parameters
        self.cds_min_length = prgArgs.cds_min_length

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
        
        # Pipeline
        self.peptide_family = prgArgs.peptide_family
        self.query_dir = prgArgs.querydir

        self.pipeline_filename = {
            '01': {'filename': '01_hmmsearch',    },
            '02': {'filename': '02_hmmsearch_seq',    },
            '03': {'filename': '03_cds',          },
            '04': {'filename': '04_signalp',  },
            '05': {'filename': '05_mature_sequences.csv', },
            '06': {'filename': '06_mature_sequences.csv', },
            '07': {'filename': '07_mature_sequences.csv', },
            '08': {'filename': '08_mature_sequences.csv', },
        }
 
        # Pipeline Properties
        self.hmm_id_dict = {}
        self.knownpep_lst = []
        self.hmm_search_files = []
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
                f.write(f">{key}\n{Fasta_Dict[key]}\n")


# --------------------------------------------------------------------------------------
def main(prgArgs):
# --------------------------------------------------------------------------------------
    
    PM_Work = PeptideMiner(prgArgs)

    #Step 01 
    read_known_peptides(PM_Work)
    hmmsearch(PM_Work)
    read_hmmsearch(PM_Work)

#==============================================================================
if __name__ == "__main__":

    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='PeptideMiner', description=" PeptideMiner")

    prgParser.add_argument("--work_dir",default=None,required=False, dest="workdir", action='store', help="Working Folder")
    prgParser.add_argument("--data_dir",default=None,required=False, dest="datadir", action='store', help="Data Folder")
    prgParser.add_argument("--cds_min_length",default=None,required=False, dest="cds_min_length", action='store', help="CDS min length")

    prgParser.add_argument("--signalp_cutoff",default=None,required=False, dest="signalp_cutoff", action='store', help="SignalP cutoff")
    prgParser.add_argument("--signalp_min_length",default=None,required=False, dest="signalp_min_length", action='store', help="SignalP min length")
    prgParser.add_argument("--signalp_path",default=None,required=False, dest="signalp_path", action='store', help="SignalP executable")

    prgParser.add_argument("--mature_min_length",default=None,required=False, dest="mature_min_length", action='store', help="Mature peptide min length")
    prgParser.add_argument("--mature_max_length",default=None,required=False, dest="mature_max_length", action='store', help="Mature peptide max length")
    prgParser.add_argument("--mature_evalue_cutoff",default=None,required=False, dest="mature_evalue_cutoff", action='store', help="Mature peptide cutoff")

    prgParser.add_argument("--peptide_family",default=None,required=True, dest="peptide_family", action='store', help="Peptide family")
    prgParser.add_argument("--query",default=None,required=True, dest="querydir", action='store', help="Query folder of fasta file")
    
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    # prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")

    prgParser.add_argument("--config",default='PeptideMiner.cfg', type=Path,is_config_file=True,help="Configuration file",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    main(prgArgs)

#==============================================================================