import os, sys
import configargparse
from pathlib import Path


import datetime
import csv

from utils.database import sql_connector
from utils.db_tasks import (upload_known_peptides, upload_hmmsearch, upload_cds, update_seqreads_signalp, 
                            upload_matureseq, upload_noduplicates, upload_annotations,
                            get_seqreads, get_noduplicates, get_summary_familyname)
from utils.hmm_tasks import run_hmmsearch,addsequence_hmmsearch,filter_hmmsearch
from utils.signalp_tasks import run_signalp
from utils.matpep_tasks import alignment, Nterm, Cterm
from utils.blast_tasks import BLASTP_QRY_HEADER, make_blastp_db, run_blastp_query, parse_blastp_query

# import pandas as pd
import numpy as np
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
        self.pipeline_csv = {
            '01': {'csv_file': '01_hmmsearch_hmm.csv',    'csv_header' : ['hmm_id','hmm_name','transcriptome']},
            '02': {'csv_file': '02_hmmsearch_seq.csv',    'csv_header' : ['seq_id','sequence','hmm_id','hmm_name']},
            '03': {'csv_file': '03_cds_seq.csv',          'csv_header' : ['cds_id','seq_id','n_cds','cds']},
            '04': {'csv_file': '04_mature_peptides.csv',  'csv_header' : ['cds_id','signalp_pos','mature_peptide']},
            '05': {'csv_file': '05_mature_sequences.csv', 'csv_header' : ['cds_id','mature_sequence']},
            '06': {'csv_file': '06_mature_sequences.csv', 'csv_header' : ['id','hmm_id','transcriptome','matseq']},
        }
 
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
        self.sql_database_file = os.path.join(self.data_dir,'sqlite.db')
        self.sql_create = os.path.join(self.data_dir,'PeptideMiner_DB.sql')
        self.db = sql_connector(data_file= self.sql_database_file, sql_create= self.sql_create)
        
    def _init_sqllite(self):
        pass

    # -------------------------------------------
    def read_known_peptides(self,Family_Name):
    # -------------------------------------------
        self.family_name = Family_Name
        self.known_peptide = []
        
        pep_dir = os.path.join(self.known_pep_dir,self.family_name)
        for k in os.listdir(pep_dir):
            if k.endswith('.fna'):
                self.known_peptide.append(os.path.join(pep_dir,k))

        for fna_file in self.known_peptide:
            if os.path.isfile(fna_file):
                logger.info(f" [Known Peptides] {self.family_name} - Reading {fna_file}")
                dict_Seq = self.read_fasta_file(fna_file)
                self.n_known_peptide = 0
                for k in dict_Seq:
                    # Parsing Fast Header
                    #>P13204_ANFB [Bos taurus]
                    #>acession_name [specie]
                    #
                    species = k.strip().split('[')[-1].replace(']','').split('(')[0]
                    accession = k.split('_')[0][1:]
                    name = ' '.join(k.split('[')[0].split('_')[1:])

                    family_id,peptide_id = upload_known_peptides(self.db, self.family_name,species,accession,name,dict_Seq[k])
                    self.n_known_peptide += 1
                    self.knownpep_lst.append({'family_id':family_id,'peptide_id':peptide_id,
                                              'species':species,'accession':accession,'name':name,
                                              'sequence':dict_Seq[k]})
        else:
            logger.error(f" [Known Peptides] {self.family_name} - {self.n_known_peptide} peptide uploaded")


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
            
            self.hmm_id_dict = {}
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

                    self.hmm_id_dict[hmm_id] = [hmm_id,hmm_name,transcriptome_name]
                    seq_id_dict[sqeq_id] = [sqeq_id,sequence,hmm_id,hmm_name]

        logger.info(f" [HMM Search] HMM's: {len(self.hmm_id_dict)} uploaded")
        logger.info(f" [HMM Search] Seq's: {len(seq_id_dict)} uploaded")
            
        # Pipeline CSV Log
        csv_filename = '01_hmmsearch_hmm.csv'
        with open(os.path.join(self.pipeline_dir,csv_filename),'w',newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csv_header=['hmm_id','hmm_name','transcriptome']
            csvwriter.writerow(csv_header)
            for key in self.hmm_id_dict:
                csvwriter.writerow(self.hmm_id_dict[key])
        logger.info(f" [HMM Search] HMM's -> {csv_filename} ({len(self.hmm_id_dict)})")

        csv_filename = '02_hmmsearch_seq.csv'
        with open(os.path.join(self.pipeline_dir,csv_filename),'w',newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csv_header=['seq_id','sequence','hmm_id','hmm_name']
            csvwriter.writerow(csv_header)
            for key in seq_id_dict:
                csvwriter.writerow(seq_id_dict[key])
        logger.info(f" [HMM Search] Seq's -> {csv_filename} ({len(seq_id_dict)})")


    # ---------------------------------------------------------
    def read_cds(self, CDS_Min_Length, Overwrite=False):

        csv_dir = self.pipeline_dir
        csv_filename = '03_cds_seq.csv'
        csv_header=['cds_id','seq_id','n_cds','cds']

        # ??? All Sequences or just from the specific hmm_id
        seq_id_dict = get_seqreads(self.db) 

        self.cds_lst = []
        for seq_id in seq_id_dict:
            
            for seq_seg in seq_id['precursor'].split('*'):
                n_cds = 0
                seq_M = seq_seg
                
                if len(seq_seg) >= CDS_Min_Length and 'M' in seq_seg:
                    seq_M = seq_seg[seq_seg.index('M'):]
                    if len(seq_M) >= CDS_Min_Length:
                        n_cds += 1                
                        self.cds_lst.append({'seq_id':seq_id['id'],'n_cds':n_cds,'cds':seq_M})
                
        for cds in self.cds_lst:
            cds_id = upload_cds(self.db,cds['cds'],cds['seq_id'])
            cds['cds_id'] = cds_id   
        logger.info(f" [HMM Search] CDS: {len(self.cds_lst)} uploded")

        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
            csvwriter.writeheader()
            for cds in self.cds_lst:
                csvwriter.writerow(cds)
        logger.info(f" [HMM Search] CDS: -> {csv_filename} ({len(self.cds_lst)} )")

        
    # ---------------------------------------------------------
    def signalp_cds(self, SignalP_Cutoff, SignalP_Min_Length,Signal_Path=None, Overwrite=False):

        csv_dir = self.pipeline_dir
        csv_filename = '04_mature_peptides.csv'
        csv_header=['cds_id','signalp_pos','mature_peptide']

        if os.path.isfile(os.path.join(csv_dir,csv_filename)):
            self.maturepep_lst = []
            with open(os.path.join(csv_dir,csv_filename)) as f:
                csvreader = csv.DictReader(f)
                for row in csvreader:
                    self.maturepep_lst.append(row)
            logger.info(f" [SignalP] MatPeptides: {csv_filename} -> ({len(self.maturepep_lst)} )")
        else:
            self.maturepep_lst = []
            if len(self.cds_lst) > 0:
                for cds in self.cds_lst:
                    signalp_pos = 0
                    mature_seq = cds['cds']
                            
                    signalp_dict = run_signalp(f"{cds['cds_id']}_{cds['n_cds']:02d}",cds['cds'],SignalP_Cutoff,Signal_Path)
    
                    if signalp_dict['SP'] == 'Y':
                        signalp_pos = int(signalp_dict['CMax_pos'])
                        mature_seq = cds['cds'][signalp_pos:]
                        if len(mature_seq) >= SignalP_Min_Length:
                            self.maturepep_lst.append({'cds_id':cds['cds_id'],'signalp_pos':signalp_pos,'mature_peptide':mature_seq})
                        else:
                            # ???? Not sure
                            self.maturepep_lst.append({'cds_id':cds['cds_id'],'signalp_pos':signalp_pos,'mature_peptide':cds['cds']})

                    update_seqreads_signalp(self.db,cds['cds_id'],signalp_pos)

            
            with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
                csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
                csvwriter.writeheader()
                for mpep in self.maturepep_lst:
                    csvwriter.writerow(mpep)
            logger.info(f" [SignalP] MatPeptides: -> {csv_filename} ({len(self.maturepep_lst)} )")


    # ---------------------------------------------------------
    def select_mature(self,E_Cutoff,Min_Length,Max_Length, Overwrite=False):
        csv_dir = self.pipeline_dir
        csv_filename = '05_mature_sequences.csv'
        csv_header=['cds_id','mature_sequence']

        self.matureseq_lst = []
        evalues = []

        for mpep in self.maturepep_lst:
            mpep_seq = mpep['mature_peptide']
            best = {'result':None, 'score': None}
            sequences = []

            for filename in self.known_peptide:
                mpep_ali = alignment(mpep_seq,filename)
                for r in mpep_ali.results:
                    evalues.append(float(r.E))

                if len(mpep_ali.results) > 0:
                    # Get best alignment and check if <E_Cutoff
                    r = mpep_ali.results[0]
                    score = r.lenseq - r.overlap
                    if float(r.E) < E_Cutoff and (best['result'] is None or best['score'] > score):
                        best['result'] = r
                        best['score']  = score
                        for r in mpep_ali.results:
                            if float(r.E) < E_Cutoff:
                                sequences.append({'before_seq':mpep_seq[:r.start_q-1],
                                                  'mid_seq':   mpep_seq[r.start_q-1:r.end_q],
                                                  'after_seq': mpep_seq[r.end_q:],
                                                  'e-value': float(r.E)})

            if len(sequences) >0:
                mature_sequences = []
                for seq in sequences:
                    nterm = Nterm(seq['before_seq'])
                    cterm = Cterm(seq['after_seq'])
                    mature_sequences.append(nterm['sequence']+seq['mid_seq']+cterm['sequence'])
                    #print(f"[{nterm['sequence']}] [{seq['mid_seq']}] [{cterm['sequence']}]")
                mature_sequences = list(set(mature_sequences))

                for i,seq in enumerate(mature_sequences):
                    if len(seq) >= Min_Length and len(seq) <= Max_Length:
                        self.matureseq_lst.append({'cds_id':mpep['cds_id'],'mature_sequence':seq})

        logger.info(f" [Fasta36] MatureSeq: E-values GeoMean: {np.exp(np.log(evalues).mean()):.5f} [Q1: {np.quantile(evalues,0.25):.5f} Mean: {np.mean(evalues):.5f} Q3: {np.quantile(evalues,0.75):.5f} ]")

        # Write CSV file
        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
            csvwriter.writeheader()
            for mpep in self.matureseq_lst:
                csvwriter.writerow(mpep)
        logger.info(f" [Fasta36] MatureSeq: CutOffs: E-value:{E_Cutoff} Length: {Min_Length}-{Max_Length}) -> {csv_filename} ({len(self.matureseq_lst)} sequences)")

    # ---------------------------------------------------------
    def upload_mature(self, Overwrite=False):

        csv_dir = self.pipeline_dir
        fasta_filename = '06_mature_sequences.fna'

        logger.info(f" [Fasta36] MatureSeq: Uploading sequences)")
        for matseq in self.matureseq_lst:
            upload_matureseq(self.db,matseq['cds_id'],matseq['mature_sequence'],verbose=1)
            upload_noduplicates(self.db,matseq['cds_id'],matseq['mature_sequence'],verbose=1)

        # Get NoDuplicates for HMM_ID
        _seq = []
        for hmm_id in self.hmm_id_dict:
            _seq += get_noduplicates(self.db,hmm_id)
        
        # Write Fasta file
        _fasta = {}
        for s in _seq:
            _name = ":".join([str(s['id']),str(s['hmm_id']),s['transcriptome']])
            _fasta[_name] = s['matseq']
        logger.info(f" [Fasta36] MatureSeq: -> {fasta_filename}")
        self.write_fasta_file(os.path.join(csv_dir,fasta_filename),_fasta)

    # ---------------------------------------------------------
    def run_blast(self, Overwrite=False):

        csv_dir = self.pipeline_dir
        csv_filename = '07_blastp_annotation.csv'
        blast_filename = '07_known_sequences'

        qry_fasta = '06_mature_sequences.fna'
        qry_out = '07_blastp_qry.csv'
        
        # Write Fasta file
        _fasta = {}
        for s in self.knownpep_lst:
            _name = ":".join([str(s['family_id']),str(s['peptide_id']),s['species'],s['name'],s['accession']])
            _fasta[_name] = s['sequence']
        logger.info(f" [BlastP] MatureSeq: -> {blast_filename}.fna")
        self.write_fasta_file(f"{os.path.join(csv_dir,blast_filename)}.fna",_fasta)

        # Make BlastP database
        make_blastp_db(os.path.join(csv_dir,blast_filename))

        # Run and Parse BlastP query for unique qry_names with lowest evalue
        run_blastp_query(os.path.join(csv_dir,blast_filename),os.path.join(csv_dir,qry_fasta),os.path.join(csv_dir,qry_out))
        blastp_out = parse_blastp_query(os.path.join(csv_dir,qry_out))
        unq_qry_names = set([b['qry_name'] for b in blastp_out])
        
        self.blastp_annotations = []
        for qn in unq_qry_names:
            _q = [i for i in blastp_out if i['qry_name'] == qn]
            _q.sort(key=lambda i:float(i['evalue']))
            self.blastp_annotations.append(_q[0])
        
        # Upload Annotation
        for a in self.blastp_annotations:
            _nid = a['qry_name'].split(':')[0]
            _kid = a['subject_name'].split(':')[0]
            upload_annotations(self.db,_nid,_kid,float(a['pct_identity']),float(a['evalue']),int(a['length']))
        logger.info(f" [BlastP] Annotations {len(self.blastp_annotations)} uploded")

        # Write CSV file
        logger.info(f" [BlastP] Annotations -> {csv_filename} ({len(self.blastp_annotations)})")
        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=BLASTP_QRY_HEADER)                
            csvwriter.writeheader()
            for a in self.blastp_annotations:
                csvwriter.writerow(a)

    # ---------------------------------------------------------
    def summary(self, FamilyName, Overwrite=False):

        csv_dir = self.pipeline_dir
        csv_filename = f'08_summary_{FamilyName}.csv'

        sum_data = get_summary_familyname(self.db,FamilyName)
        csv_header = list(sum_data[0].keys())
        # Write CSV file

        logger.info(f" [BlastP] Annotations -> {csv_filename} ({len(self.blastp_annotations)})")
        with open(os.path.join(csv_dir,csv_filename),'w',newline='') as f:
            csvwriter = csv.DictWriter(f, fieldnames=csv_header)                
            csvwriter.writeheader()
            for s in sum_data:
                csvwriter.writerow(s)




# --------------------------------------------------------------------------------------
def main(prgArgs):
# --------------------------------------------------------------------------------------
    
    pWork = PeptideMiner(WorkDir=prgArgs.workdir, 
                         DataDir=prgArgs.datadir)
    
    pWork.read_known_peptides(prgArgs.peptide_family)

    # Steps 0, 1, 2
    pWork.hmmsearch(prgArgs.querydir)
    pWork.read_hmmsearch()

    # Step 3, 4
    pWork.read_cds(int(prgArgs.cds_min_length))

    pWork.signalp_cds(float(prgArgs.signalp_cutoff),
                      int(prgArgs.signalp_min_length),
                      prgArgs.signalp_path)

    # Step 5, 6
    pWork.select_mature(float(prgArgs.mature_evalue_cutoff),
                        int(prgArgs.mature_min_length),
                        int(prgArgs.mature_max_length))
    pWork.upload_mature(prgArgs.peptide_family)

    # Step 7, 8
    pWork.run_blast()
    pWork.summary(prgArgs.peptide_family)

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

    prgParser.add_argument("--peptide_family",default=None,required=False, dest="peptide_family", action='store', help="Peptide family")
    prgParser.add_argument("--query",default=None,required=False, dest="querydir", action='store', help="Query folder of fasta file")
    
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