CREATE TABLE annotated (
  id INTEGER PRIMARY KEY,
  novel_id INTEGER DEFAULT NULL,
  knownNP_id INTEGER DEFAULT NULL,
  pid float DEFAULT NULL,
  evalue float DEFAULT NULL,
  length_alignment INTEGER DEFAULT NULL
);
CREATE TABLE cds (
  id INTEGER PRIMARY KEY,
  sequence TEXT
);
CREATE TABLE hmm (
  id INTEGER PRIMARY KEY,
  name VARCHAR DEFAULT NULL
);
CREATE TABLE known_NP (
  id INTEGER PRIMARY KEY,
  name VARCHAR DEFAULT NULL,
  familyid INTEGER DEFAULT NULL,
  species VARCHAR DEFAULT NULL,
  sequence text,
  accession VARCHAR DEFAULT NULL
);
CREATE TABLE mature (
  id INTEGER PRIMARY KEY,
  cds_id INTEGER DEFAULT NULL,
  matseq TEXT,
  noduplicates_id INTEGER DEFAULT NULL
);
CREATE TABLE neuropeptide_family (
  id INTEGER PRIMARY KEY,
  name VARCHAR DEFAULT NULL
);
CREATE TABLE noduplicates (
  id INTEGER PRIMARY KEY,
  hmm_id INTEGER DEFAULT NULL,
  transcriptome VARCHAR DEFAULT NULL,
  matseq TEXT
);
CREATE TABLE seqreads (
  id INTEGER PRIMARY KEY,
  name VARCHAR DEFAULT NULL,
  hmmid INTEGER DEFAULT NULL,
  transcriptome VARCHAR DEFAULT NULL,
  evalue float DEFAULT NULL,
  signalseq_length INTEGER DEFAULT NULL,
  precursor TEXT,
  cds_id INTEGER DEFAULT NULL
);