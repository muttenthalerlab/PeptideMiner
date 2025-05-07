import os, sys
import configargparse
from pathlib import Path


import datetime
import csv

from utils.database import sql_connector
from utils.dbtasks import upload_known_peptides, upload_hmmsearch, get_seqreads
from utils.hmmtasks import run_hmmsearch,addsequence_hmmsearch,filter_hmmsearch

# import pandas as pd
# import numpy as np
# from functools import reduce

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "PeptideMiner"
#logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
#    format="[%(name)-20s] %(message)s ",
    format="%(message)s",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)


# --------------------------------------------------------------------------------------
class PeptideMiner():
# --------------------------------------------------------------------------------------

    def __init__(self, WorkDir='.', DataDir='data'):
        
        # Data Folders
        self.data_dir = DataDir
        self.known_pep_dir = os.path.join(self.data_dir,'01-known_seq')
        self.hmm_dir = os.path.join(self.data_dir,'02-pHMM')

        # Work Folders
        self.work_dir = WorkDir
        self.hmmsearch_dir = os.path.join(self.work_dir,'01-hmmsearch')
        self.pipeline_dir = os.path.join(self.work_dir,'02-pipeline')
        
        # Pipeline
        self.hmm_search_files = []
        #self.seq_id_dict = {}

        # Initialise Working Folders
        if not os.path.exists(self.hmmsearch_dir):
            os.makedirs(self.hmmsearch_dir)
        if not os.path.exists(self.pipeline_dir):
            os.makedirs(self.pipeline_dir)

        # SQL Database
        self.sql_database_file = os.path.join(self.data_dir,'sqlite.db')
        self.sql_create = os.path.join(self.data_dir,'PeptideMiner_DB.sql')
        self.db = sql_connector(data_file= self.sql_database_file, sql_create= self.sql_create)
        
    def _init_sqllite(self):
        pass

    # -------------------------------------------
    def read_known_peptides(self,Family_Name):
    # -------------------------------------------
        self.family_name = Family_Name
        
        pep_dir = os.path.join(self.known_pep_dir,self.family_name)
        for k in os.listdir(pep_dir):
            if k.endswith('.fna'):
                FastA_File = os.path.join(pep_dir,k)

        if os.path.isfile(FastA_File):
            logger.info(f" [Known Peptides] {self.family_name} - Reading {FastA_File}")
            dict_Seq = self.read_fasta_file(FastA_File)
            self.n_known_peptide = 0
            for k in dict_Seq:
                # Parsing Fast Header
                #>P13204_ANFB [Bos taurus]
                #>acession_name [specie]
                #
                species = k.strip().split('[')[-1].replace(']','').split('(')[0]
                accession = k.split('_')[0][1:]
                name = ' '.join(k.split('[')[0].split('_')[1:])

                upload_known_peptides(self.db, self.family_name, species,accession,name,dict_Seq[k])
                self.n_known_peptide += 1
        else:
            logger.error(f" [Known Peptides] {self.family_name} - {self.n_known_peptide} peptide uploaded")


    # ----------------------------------------------------------
    @staticmethod
    def read_fasta_file(FastA_File):
        """
         Reads FastA file into Dictionary with >line as key and sequence as value
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


    # ---------------------------------------------------------
    def hmmsearch(self, Query, Overwrite=False):

        # Output Summary        
        hmm_search_outfile = os.path.join(self.hmmsearch_dir,'hmmsearchoutput.txt') 
        hmm_search_out = []
        
        # Check if hmmsearch has been run before
        self.hmm_search_files = [f for f in os.listdir(self.hmmsearch_dir) if f.endswith('csv')]
        
        if len(self.hmm_search_files)==0 or Overwrite:
            self.hmm_search_files = []
            
            #Get a list of the profile-HMMs available for the desired family
            HMM_Files = [h for h in os.listdir(self.hmm_dir) if h.startswith(self.family_name)]
            logger.info(f" [HMM Search] {self.family_name} - HMM Files : {HMM_Files}")
            
            #For each Query (fna) file run hmmsearch
            Query_Files = [q for q in os.listdir(Query) if q.endswith('fna')]
            for qry_file in Query_Files:
                
                dict_Fasta = self.read_fasta_file(os.path.join(Query,qry_file))
                logger.info(f" [HMM Search] {self.family_name} - HMM Query : {qry_file} with {len(dict_Fasta)} sequences")

                for hmm in HMM_Files:
                    hmm_out = f"{hmm.split('.')[0]}.{qry_file.replace('.fna','')}"
                    
                    # Run HMM Search 
                    _ret = run_hmmsearch(os.path.join(self.hmmsearch_dir,hmm_out),os.path.join(self.hmm_dir,hmm),os.path.join(Query,qry_file))
                    hmm_search_out.append(_ret)
                                                
                    # Add Seqeunce infor to HMM Search Output CSV file
                    addsequence_hmmsearch(self.hmmsearch_dir,hmm_out,dict_Fasta)
                    self.hmm_search_files.append(f"{hmm_out}.csv")
            
            # Output logfiles
            with open(hmm_search_outfile,'w') as out:
                for line in hmm_search_out:
                    out.writelines(line)
        else:
            logger.info(f" [HMM Search] Exists {self.hmm_search_files} in {self.hmmsearch_dir}")
                    
    # ---------------------------------------------------------
    def read_hmmsearch(self, Overwrite=False):
        
        if len(self.hmm_search_files)>0:
            logger.info(f" [HMM Search] Files: {len(self.hmm_search_files)}")
            
            hmm_id_dict = {}
            seq_id_dict = {}
                        
            for hmm in  self.hmm_search_files:
                _hmm_lst = hmm.split('.')                
                hmm_name = _hmm_lst[0]
                transcriptome_name =_hmm_lst[1]
                
                hmm_lstdict = filter_hmmsearch(self.hmmsearch_dir,hmm,hmm_name,transcriptome_name)
                for id in hmm_lstdict:
                    evalue = hmm_lstdict[id][3]
                    sequence = hmm_lstdict[id][4]
                    hmm_id, sqeq_id = upload_hmmsearch(self.db,hmm_name,transcriptome_name,id,evalue,sequence)

                    hmm_id_dict[hmm_id] = [hmm_id,hmm_name,transcriptome_name]
                    seq_id_dict[sqeq_id] = [sqeq_id,sequence,hmm_id,hmm_name]
                    
        # Pipeline CSV Log
        csv_filename = '01_hmmsearch_hmm.csv'
        with open(os.path.join(self.pipeline_dir,csv_filename),'w',newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csv_header=['hmm_id','hmm_name','transcriptome']
            csvwriter.writerow(csv_header)
            for key in hmm_id_dict:
                csvwriter.writerow(hmm_id_dict[key])
        logger.info(f" [HMM Search] HMM_ID: {len(hmm_id_dict)} -> {csv_filename}")

        csv_filename = '02_hmmsearch_seq.csv'
        with open(os.path.join(self.pipeline_dir,csv_filename),'w',newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csv_header=['seq_id','sequence','hmm_id','hmm_name']
            csvwriter.writerow(csv_header)
            for key in seq_id_dict:
                csvwriter.writerow(seq_id_dict[key])
        logger.info(f" [HMM Search] Seq_ID: {len(seq_id_dict)} -> {csv_filename}")


     # ---------------------------------------------------------
    def read_cds(self, Overwrite=False):
        seq_id_dict = get_seqreads(self.db)
        print(seq_id_dict)  

       
# --------------------------------------------------------------------------------------
def main(prgArgs):
# --------------------------------------------------------------------------------------
    
    pWork = PeptideMiner(WorkDir=prgArgs.workdir, DataDir=prgArgs.datadir)
    pWork.read_known_peptides(prgArgs.peptide_family)
    pWork.hmmsearch(prgArgs.querydir)
    pWork.read_hmmsearch()
    pWork.read_cds()


#==============================================================================
if __name__ == "__main__":

    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='PeptideMiner', description=" PeptideMiner")

    prgParser.add_argument("--work_dir",default=None,required=False, dest="workdir", action='store', help="Working Folder")
    prgParser.add_argument("--data_dir",default=None,required=False, dest="datadir", action='store', help="Data Folder")
    prgParser.add_argument("--cds_min_length",default=None,required=False, dest="cds_min_length", action='store', help="CDS min length")

    prgParser.add_argument("--sp_cutoff",default=None,required=False, dest="sp_cutoff", action='store', help="SignalP cutoff")
    prgParser.add_argument("--sp_min_length",default=None,required=False, dest="sp_min_length", action='store', help="SignalP min length")

    prgParser.add_argument("--mature_min_length",default=None,required=False, dest="mature_min_length", action='store', help="Mature peptide min length")
    prgParser.add_argument("--mature_max_length",default=None,required=False, dest="mature_max_length", action='store', help="Mature peptide max length")
    prgParser.add_argument("--mature_evalue_cutoff",default=None,required=False, dest="mature_evalue_cutoff", action='store', help="Mature peptide cutoff")

    prgParser.add_argument("--peptide_family",default=None,required=False, dest="peptide_family", action='store', help="Peptide family")
    prgParser.add_argument("--query",default=None,required=False, dest="querydir", action='store', help="Query folder of fasta file")
    
    # prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [TestPlate]")
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    # prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    # prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--test",default=0,required=False, dest="te
    # st", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    #prgParser.add_argument("-p","--pivot",default=None,required=False, dest="pivot", action='store', help="Pivot tables")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    # prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
    # prgParser.add_argument("-e","--excel",default=None,required=False, dest="excelfile", action='store', help="Excel File")
    # prgParser.add_argument("-f","--format",default='Check',required=False, dest="pivot", action='store', help="Format of output EXcel")
    # prgParser.add_argument("--plotdir",default=None,required=False, dest="plotdir", action='store', help="Folder for Plots")
    # #prgParser.add_argument("-o","--outdir",default=None,required=False, dest="outdir", action='store', help="Prefix to add to PlateID")

    prgParser.add_argument("-c","--config",default='PeptideMiner.cfg', type=Path,is_config_file=True,help="Configuration file",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    main(prgArgs)

#==============================================================================