import re,sys,os,random
import subprocess

"""
Step 4 support module.
Module to connect with and run SignalP for step 4.
Currently accomodates SignalP 4.1 (http://www.cbs.dtu.dk/services/SignalP-4.1/)
"""

def create_temp_fasta (sequence,name='seq'):
	tmp_seq_file = f"/tmp/seq_{random.randint(1,10e9)}.fasta"
	tmp = open(tmp_seq_file,'w')
	tmp.write(f">{name}\n{sequence}")
	tmp.close()
	return tmp_seq_file

def find (sequence,signalp,cutoff):
    #Run signalp; for signalp-4.1"""
	tmp_seq_file = create_temp_fasta(sequence)
    # cmd = f"{signalp} -t euk -U {cutoff} -u {cutoff} {tmp_seq_file}".format(signalp,,cutoff)"
    # p = subprocess.run(cmd,shell=True,capture_output=False, text=True)

	signalp = os.popen("{0} -t euk -U {2} -u {2} {1}".format(signalp,tmp_seq_file,cutoff))
	
	#Parse output"""
	pos = 0
	insignal = 0
	inHMM = 0
	linepat = re.compile("^[^ ]+ +[0-9\\.]+ +(\\d+) +[0-9\\.]+ +\\d+ +[0-9\\.]+ +\\d+ +[0-9\\.]+ +[0-9\\.]+ +([YN]) +[0-9\\.]+")
	for line in signalp:
		if line[0] == "#": continue
		s = linepat.search(line)
		if s:
			if s.group(2) == "Y":
				pos = int(s.group(1))-1
			break
	signalp.close()
	os.system("rm {0}".format(tmp_seq_file))
	return pos


