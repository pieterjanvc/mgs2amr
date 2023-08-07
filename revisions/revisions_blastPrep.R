#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)
library(gfaTools)

args = commandArgs(trailingOnly = TRUE)
baseFolder = args[[1]] %>% gfaTools::formatPath()
pipelineId = args[[2]] %>% as.integer()
resultsDB = args[[3]]

trim = 5000 #max environment around an ARG

overlap = function(x, y, amount = 50){
  grp = 1
  if(length(x) == 1) return(grp)
  out = rep(1, length(x))
  z = x
  x = ifelse(x > y, y, x)
  y = ifelse(z < y, y, z)
  x = c(x, x[length(x)])
  y = c(y, y[length(y)])
  for(i in 1:(length(x) - 2)){
    if(min(y[i],y[i+1]) - max(x[i], x[i+1]) >= amount){
      out[i+1] = grp
    } else {
      grp = grp + 1
      out[i+1] = grp
    }
  }
  return(out)
}

#Get pipeline info from the database
myConn = dbConnect(dbDriver("SQLite"), resultsDB)
pipelineInfo = tbl(myConn, "pipeline") %>% filter(pipelineId == {{pipelineId}}) %>% collect()
dbDisconnect(myConn)

#Load data from diamond
diaData = read.csv(paste0(pipelineInfo$folder, "/diamondOutput.csv"), sep = "\t", header = F) 
colnames(diaData) = c("qseqid", "sseqid", "qlen", "slen", "pident", "nident",
                      "length", "evalue", "bitscore", "qstart", "qend")

myConn = dbConnect(dbDriver("SQLite"), paste0(baseFolder, "SupplementalDatabase.db"))
ARG = tbl(myConn, "ARG") %>% collect()
bactInfo = tbl(myConn, "genotypes") %>% filter(Run == local(pipelineInfo$isolate)) %>% 
  left_join(tbl(myConn, "geneLinks"), by = "gene") %>% collect() %>% distinct() %>% 
  filter(!is.na(geneId)) %>% group_by(geneId) %>% 
  filter(type == ifelse("COMPLETE" %in% type, "COMPLETE", type[1])) %>% ungroup()
dbDisconnect(myConn)

argInfo = diaData %>% 
  left_join(ARG %>% select(geneId, prot, gene, subtype), by = c("sseqid" = "prot")) %>% 
  mutate(coverage = nident / slen) %>% filter(coverage >= 0.5) %>%
  left_join(bactInfo %>% select(geneId, type), by = "geneId") %>% 
  filter(!is.na(gene)) %>% 
  group_by(qseqid, qstart, qend, gene) %>% filter(evalue == min(evalue)) %>% 
  group_by(qseqid, gene) %>% arrange(qstart, qend) %>% 
  mutate(grp = overlap(qstart, qend)) %>% 
  group_by(qseqid, gene, grp) %>% filter(bitscore == max(bitscore)) %>% ungroup()

#Add to results database
myConn = dbConnect(dbDriver("SQLite"), resultsDB)
q = dbSendStatement(myConn, "DELETE FROM argInfo WHERE pipelineId == ?",
                    params = pipelineId)
q = dbClearResult(q)

q = dbSendStatement(myConn, 
  "INSERT INTO argInfo (qseqid, geneId, sseqid, qlen, slen, pident, nident, length, 
  evalue, bitscore, qstart, qend, gene, subtype, coverage, type, grp, pipelineId)
  VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?)", 
  params = argInfo %>% 
    select(qseqid, geneId, sseqid, qlen, slen, pident, nident, length,
           evalue, bitscore, qstart, qend, gene, subtype, coverage, type, grp) %>% 
    mutate(pipelineId = {{pipelineId}}) %>%
    as.list() %>% unname())
q = dbClearResult(q)

dbDisconnect(myConn)


#Get the scaffolds, filter and crop them
scaffolds = fasta_read(paste0(pipelineInfo$folder, "/scaffolds.fasta"), "n") %>% 
  filter(id %in% argInfo$qseqid) %>% 
  extract(id, c("node", "length", "coverage"), 
          "(NODE_\\d+)_length_(\\d+)_cov_(\\d+.?\\d+)", remove = F) %>% 
  mutate(length = as.integer(length), coverage = as.numeric(coverage) %>% round(2)) %>% 
  left_join(argInfo %>% select(qseqid, qstart, qend, geneId), by = c("id" = "qseqid")) 

scaffolds[scaffolds$qstart > scaffolds$qend, c("qstart", "qend")] = 
  scaffolds[scaffolds$qstart > scaffolds$qend, c("qend", "qstart")]
  
#Add to results database
myConn = dbConnect(dbDriver("SQLite"), resultsDB)
q = dbSendStatement(myConn, "DELETE FROM scaffold WHERE pipelineId == ?",
                    params = pipelineId)
q = dbClearResult(q)

q = dbSendStatement(myConn, 
  "INSERT INTO scaffold (scaffoldId, geneId, length, coverage, geneStart, geneEnd, pipelineId)
  VALUES(?, ?, ?, ?, ?, ?, ?)", 
  params = scaffolds %>% 
    select(id, geneId, length, coverage, qstart, qend) %>% 
    mutate(pipelineId = {{pipelineId}}) %>%
    as.list() %>% unname())
q = dbClearResult(q)
dbDisconnect(myConn)

#save for BLASTn as new fasta
scaffolds = scaffolds %>% rowwise() %>% 
  mutate(
    geneSeq = str_sub(seq, max(qstart - trim, 1), min(qend + trim, length)),
    id = paste(id, geneId, sep = "_")
  )

fasta_write(scaffolds$geneSeq, paste0(pipelineInfo$folder, "/toBlast.fasta"), 
            id = scaffolds$id, checkData = F)
