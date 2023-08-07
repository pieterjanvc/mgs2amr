#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)

args = commandArgs(trailingOnly = TRUE)
pipelineId = args[[1]] %>% as.integer()
resultsDB = args[[2]]

myConn = dbConnect(dbDriver("SQLite"), resultsDB)
pipelineInfo = tbl(myConn, "pipeline") %>% filter(pipelineId == {{pipelineId}}) %>% collect()
dbDisconnect(myConn)

#Read the BLASTn results
results = read_delim(paste0(pipelineInfo$folder, "/blastResults.csv"), 
			delim = "\t", col_names = F, show_col_types = FALSE)
colnames(results) = c("qseqid", "sallacc", "staxids", "sscinames", "salltitles",
                      "qlen", "slen", "qstart", "qend", "sstart", "send",
                      "bitscore", "score", "length", "pident", "nident",
                      'qcovs', 'qcovhsp')

#Filter to keep only the best per species / scaffold
result = results %>% group_by(qseqid, sscinames) %>% 
  filter(bitscore == max(bitscore)) %>% slice(1) %>% 
  group_by(qseqid) %>% mutate(
    rank = rank(1/bitscore, ties.method = "min"),
    singleTop = n_distinct(sscinames) == 1) %>% 
  arrange(qseqid, rank) %>% ungroup()

#Save to database
myConn = dbConnect(dbDriver("SQLite"), resultsDB)
q = dbSendStatement(myConn, "DELETE FROM blastOut WHERE pipelineId == ?",
                    params = pipelineId)
q = dbClearResult(q)

q = dbSendStatement(myConn, 
  "INSERT INTO blastOut (qseqid, sallacc, staxids, sscinames, salltitles, qlen, slen, qstart, 
  qend, sstart, send, bitscore, length, pident, nident, qcovs, rank, singleTop, pipelineId)
  VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", 
  params = result %>% select(
	qseqid, sallacc, staxids, sscinames, salltitles, qlen, slen, qstart, qend, 
	sstart, send, bitscore, length, pident, nident, qcovs, rank, singleTop) %>% 
	mutate(pipelineId = {{pipelineId}}) %>% as.list() %>% unname())
q = dbClearResult(q)
dbDisconnect(myConn)
