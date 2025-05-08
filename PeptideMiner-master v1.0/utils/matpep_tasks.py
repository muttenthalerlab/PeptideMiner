import re,glob,random,gzip,bz2,os

import os
import csv
import subprocess

import logging
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------
def find_cys_region (sequence):
# -----------------------------------------------------------------------
    pre_seq = ''
    post_seq = ''
    m = re.search(r"^([^C]*)(C.*C)(.*)",sequence)
    if m:
        pre_seq = m.group(1)
        sequence = m.group(2)
        post_seq = m.group(3)
    return pre_seq,sequence,post_seq

# -----------------------------------------------------------------------
def find_regions (sequence,fasta,dir_fasta,ecutoff):
# -----------------------------------------------------------------------
    best_r = None
    best_score = None
    for filename in glob.glob("{0}/*.fna*".format(dir_fasta)):
        name = re.search(r"/(.*)\.fna",filename).group(1)
        ali = alignment(sequence,fasta,filename)
        if len(ali.results) == 0: 
            continue
        r = ali.results[0]
        score = r.lenseq - r.overlap
        if float(r.E) < ecutoff and (best_r is None or best_score > score):
            best_r = r
            best_score = score
            sequences = []
            for r in ali.results:
                if float(r.E) < float(ecutoff):
                    sequences.append([sequence[:r.start_q-1],sequence[r.start_q-1:r.end_q],sequence[r.end_q:]])
    if best_r is None: return []
    return sequences

def Nterm (sequence,use_isolated=1):
	cleavages = []
	Furin_cleav  = re.compile("^(.*)([KR][KR].R)(.*)$")
	RRGC_cleav   = re.compile("^(.*)(RR.)(G)$")
	KRR_cleav    = re.compile("^(.*)(KR)(R.*)$")
	LVLK_cleav   = re.compile("^(.*)([LV]..L..K)(.*)$")
	LSR_cleav    = re.compile("^(.*)(L..S.R)(.*)$")
	LENDKR_cleav = re.compile("^(.*)(L.[END]KR)(.*)$")
	LKR_cleav    = re.compile("^(.*)(L..[KR])(.*)$")
	KR_cleav     = re.compile("^(.*)(KR)(.*)$")
	RR_cleav     = re.compile("^(.*)(RR)(.*)$")
	KK_cleav     = re.compile("^(.*)(KK)(.*)$")
	KE_cleav     = re.compile("^(.*)(KE)(.?)$")
	ER_cleav     = re.compile("^(.*)(ER)(.*)$")
	E_R_cleav    = re.compile("^(.*)(E...R)(.*)$")
	E__R_cleav   = re.compile("^(.*)(E.....R)(.*)$")
	RKtoR_cleav  = re.compile("^(.*)([RK].{0,5}R)(.*)$")
	DtoR_cleav   = re.compile("^(.*)(D.{0,6}R)(.*)$")
	RMVL_cleav   = re.compile("^(.*)([RM].VL)(.*)$")

	m = Furin_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('Furin',m))
		sequence = m.group(3)
	m = RRGC_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('RR.GC',m))
		sequence = 'G'
	# GIIIA impose un cleavage KR|R et donc montre une plus grande affinite sur KR
	# que sur RR (ou alors deux enzymes et une cynetiquement plus rapide).
#	m = KRR_cleav.search(sequence) # new
#	if m:
#		cleavages.append(cleavage('KRR',m))
#		sequence = m.group(3)
#	m = KR_cleav.search(sequence) # new
#	if m:
#		cleavages.append(cleavage('KR',m))
#		sequence = m.group(3)

	# PIIIA proves that the lowest number of residues between RK and R is selected
	# GIIIA proves that K is prefered to R

	m = LVLK_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('[LV]..L..K',m))
		sequence = m.group(3)
	m = LSR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('L..S.R',m))
		sequence = m.group(3)
	m = LENDKR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('L.[END]KR',m))
		sequence = m.group(3)
	m = LKR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('L..[KR]',m))
		sequence = m.group(3)
	m = KR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('KR',m))
		sequence = m.group(3)
	m = RR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('RR',m))
		sequence = m.group(3)
	m = KK_cleav.search(sequence)
	if m and use_isolated:
		cleavages.append(cleavage('KK',m))
		sequence = m.group(3)
	m = KE_cleav.search(sequence)
	if m and use_isolated:
		cleavages.append(cleavage('KE',m))
		sequence = m.group(3)

	# EbeforeR find with spacing 3 and 5
	m = ER_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('ER',m))
		sequence = m.group(3)
	m = E_R_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('E...R',m))
		sequence = m.group(3)
	m = E__R_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('E.....R',m))
		sequence = m.group(3)

	m = RKtoR_cleav.search(sequence)
	if m:
		if m.group(2) == 'RPR':
			pre = "{0}{1}".format(m.group(2),m.group(3))
			tp = m.group(1)
			m2 = RKtoR_cleav.search(tp)
			if m2: # remove RPR case for Pl14a
				cleavages.append(cleavage('R(K)toR',m))
				sequence = "{0}{1}".format(m2.group(3),pre)
		else:
			cleavages.append(cleavage('R(K)toR',m))
			sequence = m.group(3)

	m = DtoR_cleav.search(sequence)
	if m: # new one only for MrIIIG
		cleavages.append(cleavage('DtoR',m))
		sequence = m.group(3)

	m = RMVL_cleav.search(sequence)
	if m and use_isolated:
		cleavages.append(cleavage('[RM].VL',m))
		sequence = m.group(3)

	return {'sequence':sequence,'cleavages':cleavages}

def Cterm (sequence,use_isolated=1):
	cleavages = []
	KR_cleav    = re.compile("(.*?)([KR][KR])(.*)")
	RTIL_cleav  = re.compile("(.*?)(RT[IL])(.*)$")
#	ER_cleav    = re.compile("^(.*)ER(.*)$")
#	RKtoR_cleav = re.compile("^(.*)([RK].{0,6}R)(.*)$")
#	RKtoK_cleav = re.compile("^(.*)([RK].{0,6}K)(.*)$")
#	EtoR_cleav  = re.compile("^(.*)(E.{0,6}R)(.*)$")
#	RtoE_cleav  = re.compile("^(.*)([RK].{0,6}E)(.*)$")

	m = KR_cleav.search(sequence)
	if m:
		cleavages.append(cleavage('KR',m))
		sequence = m.group(1)
	m = RTIL_cleav.search(sequence)
	if m and use_isolated:
		cleavages.append(cleavage('RT[IL]',m))
		sequence = m.group(1)
# Commenting those lines or MIVA cleavage
#	m = ER_cleav.search(sequence)
#	if m:
#		cleavages.append(cleavage('ER',m))
#		sequence = m.group(1)
#	m = RKtoR_cleav.search(sequence)
#	if m:
#		cleavages.append(cleavage('R(K)toR',m))
#		sequence = m.group(1)
#	m = RKtoK_cleav.search(sequence)
#	if m:
#		cleavages.append(cleavage('R(K)toK',m))
#		sequence = m.group(1)
#	m = EtoR_cleav.search(sequence)
#	if m:
#		cleavages.append(cleavage('EbeforeR',m))
#		sequence = m.group(1)
#	m = RtoE_cleav.search(sequence)
#	if m:
#		cleavages.append(cleavage('EafterR',m))
#		sequence = m.group(1)
	return {'sequence':sequence,'cleavages':cleavages}

def CPE (sequence):
	CPE_cleav   = re.compile("(.+[^KR])([KR]+)$")
	m = CPE_cleav.search(sequence)
	if m:
		return m

def PAM (sequence):
	PAM_cleav   = re.compile("(.+)(G)$")
	m = PAM_cleav.search(sequence)
	if m:
		return m

class cleavage:
	def __init__ (self,name,m):
		self.name = name
		self.sequence = m.group(2)
		self.pre = m.group(1)
		self.post = m.group(3)

	def __str__ (self):
		return(f"  {self.name} cleavage: {self.pre} - {self.sequence} - {self.post}")

def create_temp_fasta (sequence,name='seq'):
	tmp_seq_file = f"/tmp/seq_{random.randint(1,10e9)}.fasta"
	tmp = open(tmp_seq_file,'w')
	tmp.write(f">{name}\n{sequence}")
	tmp.close()
	return tmp_seq_file

# -----------------------------------------------------------------------
class alignment:
# -----------------------------------------------------------------------
    def __init__ (self,qry_sequence,lib_file='',gap_penalty=-10,evalue='1'):
        self.gap_penalty = gap_penalty
        self.evalue = evalue
        self.qry_sequence = qry_sequence
        self.lib_file = lib_file
		
        self.fasta36_out = []
        self.results = []

        self.run_fasta36()
        self.read_results()
        
    def run_fasta36(self):
        tmp_seq_file = create_temp_fasta(self.qry_sequence)
        cmd = f"fasta36 -f {self.gap_penalty} -z 0 -E {self.evalue} -q {tmp_seq_file} {self.lib_file} "
        #print(cmd)
        p = subprocess.run(cmd,shell=True,capture_output=True, text=True, encoding='latin-1')
        self.fasta36_out = p.stdout.splitlines()
        # print(self.fasta36_out)
        os.remove(tmp_seq_file)
		
    def __str__ (self):
        res_str = f"{len(self.results)} Result(s)"
        for i in range(0,len(self.results)):
            res_str += f"\n [{i:02d}] : {self.results[i].summary()}"
        return res_str

    def __iter__(self):
        return self.results.__iter__()

    def read_results (self):
        record = 0
        res = None
        for line in self.fasta36_out:
            line=line.rstrip()
            if re.search(r"\d+ residues in",line):
                record = 0
            if re.search(r"^>>",line):
                if res is not None:
                    res.parse_alignment()
                res = result(line.lstrip(">"))
                self.results.append(res)
                record = 1
            elif re.search(r"^>--",line):
                if res is not None:
                    res.parse_alignment()
                res = result(name=res.name,dbid=res.dbid,lenseq=res.lenseq)
                self.results.append(res)
                record = 1
            elif record and res.score1 == None:
                res.score1 = line
                res.parse_score1(line)
            elif record and res.score2 == None:
                res.score2 = line
                res.parse_score2(line)
            elif record:
                res.ali = f"{res.ali}\n{line}"
        if res is not None:
            res.parse_alignment()

# -----------------------------------------------------------------------
class result:
# -----------------------------------------------------------------------
    def __init__ (self,linename=None,name=None,dbid=None,lenseq=0):
        if linename is not None:
            m = re.search(r'(.*[^ ]) +\((\d+) (aa|nt)\)',linename)
            name = m.group(1)
            lenseq = int(m.group(2))
            m = re.search(r'^(\d+) (.*)',name)
            if m:
                self.dbid    = int(m.group(1))
                self.name    = m.group(2)
            else:
                self.dbid    = None
                self.name    = name
        elif name is not None:
            self.dbid    = dbid
            self.name    = name
        self.lenseq  = lenseq
        self.score1  = None
        self.score2  = None
        self.Z       = None
        self.E       = 57
        self.SW      = None
        self.PID     = None
        self.PS      = None
        self.overlap = None
        self.start_q = None
        self.end_q   = None
        self.start_l = None
        self.end_l   = None
        self.sense   = None
        self.frame   = None
        self.ali     = ''
        self.ali1    = ''
        self.ali2    = ''
        self.alim    = ''

    def parse_score1 (self,line):
        m = re.search(r"initn: +([0-9]+) init1: +([0-9]+) opt: +([0-9]+) +Z-score: +([0-9\.]+) +bits: +([0-9\.]+) +E\([^\)]*\): +([0-9\.eE\-\+]+)",line)
        if m:
            self.Z = m.group(4)
            self.E = m.group(6)
        m = re.search(r"Frame: (.) initn",line)
        if m:
            self.sense = m.group(1)

    def parse_score2 (self,line):
        m = re.search(r"Smith-Waterman score: (\d+); ([0-9\.]+)% identity \(([0-9\.]+)% similar\) in (\d+) (aa|nt) overlap \((\d+)-(\d+):(\d+)-(\d+)\)",line)
        if m:
            self.SW      = float(m.group(1))
            self.PID     = float(m.group(2))
            self.PS      = float( m.group(3))
            self.overlap = int( m.group(4))
            self.start_q = int( m.group(6))
            self.end_q   = int( m.group(7))
            self.start_l = int( m.group(8))
            self.end_l   = int( m.group(9))
            if self.sense is not None:
                if self.sense == 'f':
                    self.frame = ((self.start_l - 1)% 3) + 1
                else:
                    self.frame = -1*((self.lenseq-self.start_l) % 3 + 1)

    def parse_alignment(self):
        lines = self.ali.split("\n")
        ali1 = []
        ali2 = []
        alim = []
        if len(lines) > 3:
            nb_clust_ali = int((len(lines)-2)/6)
            for i in range(0,nb_clust_ali):
                linenum = i*6 + 2
                a1 = lines[linenum+1][7:]
                am = lines[linenum+2][7:]
                a2 = lines[linenum+3][7:]
                if i == nb_clust_ali - 1:
                    a1 = a1.lstrip()
                    a2 = a2.lstrip()
                    minlen = min(len(a1),len(a2))
                    a1 = a1[0:minlen]
                    am = am[0:minlen]
                    a2 = a2[0:minlen]
                else:
                    am = "{0}{1}".format(am," "*(len(a1)-len(am)))
                if i == 0:
                    blanksearch = re.search("^( *)",am)
                    nbblank = len(blanksearch.group(1))
                    a1 = a1[nbblank:]
                    am = am[nbblank:]
                    a2 = a2[nbblank:]
                ali1.append(a1)
                ali2.append(a2)
                alim.append(am)
        self.ali1 = "".join(ali1)
        self.ali2 = "".join(ali2)
        self.alim = "".join(alim)

    def summary(self):
        return(f"[{self.PID:3.0f} PID - {self.PS:3.0f} PS] [{self.overlap:3d} Overlap] {self.name}")

    def __str__ (self):
        retStr  = f">{self.name}\n"
        retStr += f"{self.score1}\n"
        retStr += f"{self.score2}\n"
        retStr += f"SW {self.SW} - PID {self.PID} - PS {self.PS} "
        retStr += f"Overlap {self.overlap} - Query [{self.start_q} - {self.end_q}] "
        retStr += f"| Subject [{self.start_l} - {self.end_l}] \n"
        retStr += f"E {self.E} Z {self.Z}\n"
        retStr += f"{self.ali}"
        return(retStr)