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

#from lib.database import sql_connector
from lib.pipeline_tasks import PeptideMiner, read_known_peptides, read_cds, summary
from lib.hmm_tasks import hmmsearch, read_hmmsearch
from lib.signalp_tasks import signalp_cds
from lib.matpep_tasks import select_mature, upload_mature
from lib.blast_tasks import run_blast

# --------------------------------------------------------------------------------------
def main(prgArgs):
# --------------------------------------------------------------------------------------
    
    PM_Work = PeptideMiner(prgArgs)

    err_Dict = PM_Work.check_paths()
    if err_Dict:
        print(err_Dict)
    else:
        logger.info(f" [PeptideMiner] Peptide Family: {PM_Work.family_name} - Query {PM_Work.query_dir}")
        #Step 1, 2 
        read_known_peptides(PM_Work)
        hmmsearch(PM_Work)
        read_hmmsearch(PM_Work)

        #Step 3,4
        read_cds(PM_Work)
        signalp_cds(PM_Work)


        # Step 5, 6
        select_mature(PM_Work)
        upload_mature(PM_Work)

        # Step 7, 8
        run_blast(PM_Work)
        summary(PM_Work)


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