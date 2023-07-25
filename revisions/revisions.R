#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

library(tidyverse)
library(RSQLite)
library(gfaTools)

args = commandArgs(trailingOnly = TRUE)
baseFolder = args[[1]] %>% gfaTools::formatPath()
isolate = args[[2]]
background = args[[3]]
ra = args [[4]]


isolate="SRR3345827"
background="ERR3277283"
projectDB = "D:/GoogleDrive/My Drive/Varia/UC/SupplementalDatabase.db"
baseFolder=  "D:/Desktop/"
trim = 5000 #max environment around an ARG
ra = "0.028"
msFolder = sprintf("%s/MSout/%s_%s_%s/", baseFolder, isolate, background, ra)


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

diaData = read.csv(paste0(baseFolder, "diamond/diamondOutput.csv"), sep = "\t", header = F) 
colnames(diaData) = c("qseqid", "sseqid", "qlen", "slen", "pident", "nident",
                      "length", "evalue", "bitscore", "qstart", "qend")

myConn = dbConnect(dbDriver("SQLite"), paste0(baseFolder, "SupplementalDatabase.db"))
ARG = tbl(myConn, "ARG") %>% collect()
bactInfo = tbl(myConn, "genotypes") %>% filter(Run == {{isolate}}) %>% 
  left_join(tbl(myConn, "geneLinks"), by = "gene") %>% collect() %>% 
  filter(!is.na(geneId))
bactName = tbl(myConn, "sampleInfo") %>% filter(Run == {{isolate}}) %>% pull(Bacterium)
dbDisconnect(myConn)

argInfo = diaData %>% 
  left_join(ARG %>% select(geneId, prot, gene, subtype), by = c("sseqid" = "prot")) %>% 
  mutate(coverage = round(nident / slen, 3)) %>% filter(coverage >= 0.5) %>%
  left_join(bactInfo %>% select(geneId, type), by = "geneId") %>% 
  filter(!is.na(gene)) %>% 
  group_by(qseqid, qstart, qend, gene) %>% filter(evalue == min(evalue)) %>% 
  group_by(qseqid, gene) %>% arrange(qstart, qend) %>% 
  mutate(grp = overlap(qstart, qend)) %>% 
  group_by(qseqid, gene, grp) %>% filter(bitscore == max(bitscore)) %>% ungroup()

# argInfo %>% 
#   select(qseqid, geneId, sseqid, qlen, slen, pident, nident, length,
#          evalue, bitscore, qstart, qend, gene, subtype, coverage, type, grp) %>% 
#   mutate(pipelineId = paste(isolate, background, ra, sep = "_"))

myConn = dbConnect(dbDriver("SQLite"), paste0(baseFolder, "results.db"))
# dbListFields(myConn, "argInfo") %>% paste(collapse = ", ")
q = dbSendStatement(myConn, "DELETE FROM argInfo WHERE pipelineId == ?",
                    params = paste(isolate, background, ra, sep = "_"))
q = dbClearResult(q)

q = dbSendStatement(myConn, 
  "INSERT INTO argInfo (qseqid, geneId, sseqid, qlen, slen, pident, nident, length, 
  evalue, bitscore, qstart, qend, gene, subtype, coverage, type, grp, pipelineId)
  VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?)", 
  params = argInfo %>% 
    select(qseqid, geneId, sseqid, qlen, slen, pident, nident, length,
           evalue, bitscore, qstart, qend, gene, subtype, coverage, type, grp) %>% 
    mutate(pipelineId = paste(isolate, background, ra, sep = "_")) %>%
    as.list() %>% unname())
q = dbClearResult(q)
dbDisconnect(myConn)

#Get the scaffolds, filter and cropt them
scaffolds = fasta_read(paste0(msFolder, "scaffolds.fasta"), "n") %>% 
  filter(id %in% argInfo$qseqid) %>% 
  extract(id, c("node", "length", "coverage"), 
          "(NODE_\\d+)_length_(\\d+)_cov_(\\d+.?\\d+)", remove = F) %>% 
  mutate(length = as.integer(length), coverage = as.numeric(coverage) %>% round(2)) %>% 
  left_join(argInfo %>% select(qseqid, qstart, qend, geneId), by = c("id" = "qseqid")) 

scaffolds[scaffolds$qstart > scaffolds$qend, c("qstart", "qend")] = 
  scaffolds[scaffolds$qstart > scaffolds$qend, c("qend", "qstart")]

myConn = dbConnect(dbDriver("SQLite"), paste0(baseFolder, "results.db"))
# dbListFields(myConn, "blastOut") %>% paste(collapse = ", ")
q = dbSendStatement(myConn, "DELETE FROM scaffold WHERE pipelineId == ?",
                    params = paste(isolate, background, ra, sep = "_"))
q = dbClearResult(q)

q = dbSendStatement(myConn, 
  "INSERT INTO scaffold (scaffoldId, geneId, length, coverage, geneStart, geneEnd, pipelineId)
  VALUES(?, ?, ?, ?, ?, ?, ?)", 
  params = scaffolds %>% 
    select(id, geneId, length, coverage, qstart, qend) %>% 
    mutate(pipelineId = paste(isolate, background, ra, sep = "_")) %>%
    as.list() %>% unname())
q = dbClearResult(q)
dbDisconnect(myConn)

scaffolds = scaffolds %>% rowwise() %>% 
  mutate(
    geneSeq = str_sub(seq, max(qstart - trim, 1), min(qend + trim, length)),
    id = paste(id, geneId, sep = "_")
  )

#save for BLASTn as new fasta
# fasta_write(scaffolds$geneSeq, paste0(msFolder, "toBlast.fasta"), 
#             id = scaffolds$id, checkData = F)

#CHeck after blastn
utils::browseURL(sprintf("https://www.ncbi.nlm.nih.gov/pathogens/isolates#%s", 
                         isolate))

results = read_delim("D:/Desktop/blastResults.csv", delim = "\t", col_names = F,
                     show_col_types = FALSE)
colnames(results) = c("qseqid", "sallacc", "staxids", "sscinames", "salltitles",
                      "qlen", "slen", "qstart", "qend", "sstart", "send",
                      "bitscore", "score", "length", "pident", "nident",
                      'qcovs', 'qcovhsp')

result = results %>% group_by(qseqid, sscinames) %>% 
  filter(bitscore == max(bitscore)) %>% slice(1) %>% 
  group_by(qseqid) %>% mutate(
    rank = rank(1/bitscore, ties.method = "min"),
    singleTop = n_distinct(sscinames) == 1) %>% 
  arrange(qseqid, rank) %>% ungroup()

result = result %>% filter(sscinames %in% {{bactName}}) %>% 
  group_by(qseqid) %>% slice(1) %>% ungroup() %>% 
  extract(qseqid, c("qseqid", "geneId"), "^(.*)_(\\d+$)") %>% 
  mutate(geneId = as.integer(geneId))

argInfo = argInfo %>% left_join(result %>% select(qseqid, geneId, rank, singleTop), 
                                by = c("qseqid", "geneId")) 

result %>% select(
  qseqid, sallacc, staxids, sscinames, salltitles, qlen, 
  slen, qstart, qend, sstart, send, bitscore, length, pident, nident, qcovs, rank, singleTop
) %>% mutate(pipelineId = {{pipelineId}})

#B864AT1601N


#---- old code
#--------------

#Read in the full GFA
gfa = gfa_read("D:/Desktop/assembly_graph_with_scaffolds.gfa")

#Link the qsedId (scaffold) to the node paths in the assembly graph
allPaths = gfa$paths %>% 
  mutate(pathName  = str_remove(pathName, "_\\d+$")) %>% 
  filter(pathName %in% argInfo$qseqid) %>% 
  group_by(pathName) %>% 
  summarise(segmentNames = paste(segmentNames, collapse = ","), .groups = "drop")

allPaths = setNames(lapply(allPaths$segmentNames, function(x){
  str_extract_all(x, pattern = "\\d+")[[1]]}), allPaths$pathName)

gfa = gfa[1:2] #path info in gfa no longer needed now 

# names(allPaths) = nodeNames[x]



#Get the membership of each component (subgraph id)
membership = igraph::graph_from_data_frame(data.frame(
  from = gfa$links$from,
  to = gfa$links$to
), directed = F)

membership = igraph::components(membership)$membership
membership = data.frame(
  name = names(membership),
  group = membership,
  row.names = NULL
)

#Only keep fragments which belong to relevant subgraphs (containing  GFA)
membership = membership %>% 
  filter(group %in% (membership %>% filter(name %in% sapply(allPaths, "[[", 1)) %>% 
                       pull(group))) %>% pull(name)
gfa = gfa_filterSegments(gfa, c(membership, allPaths %>% unlist()))
gfa$segments = gfa$segments %>% mutate(
  LN = nchar(sequence)
)



#For each gene, create a separate GFA to match SEQ2MGS - overlap = 55
# "NODE_431_length_11029_cov_261.585293" is loops issue

x = 3
x = which(names(allPaths) == "NODE_2096_length_2708_cov_22.146250")

myARG = allPaths[[x]] %>% unique()
names(allPaths)[x]
myARG %>% paste(collapse = ",")
#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

library(tidyverse)
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
bactInfo = tbl(myConn, "genotypes") %>% filter(Run == {{pipelineInfo$isolate}}) %>% 
  left_join(tbl(myConn, "geneLinks"), by = "gene") %>% collect() %>% 
  filter(!is.na(geneId))
dbDisconnect(myConn)

argInfo = diaData %>% 
  left_join(ARG %>% select(geneId, prot, gene, subtype), by = c("sseqid" = "prot")) %>% 
  mutate(coverage = nident / slen) %>% filter(coverage >= 0.5) %>%
  left_join(bactInfo %>% select(geneId, type), by = "geneId") %>% 
  filter(!is.na(gene)) %>% 
  group_by(qseqid, qstart, qend, gene) %>% filter(evalue == min(evalue)) %>% 
  group_by(qseqid) %>% arrange(qstart, qend) %>% 
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
                    params = paste(isolate, background, ra, sep = "_"))
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

newGfa = gfa_updateNames(gfa, myARG, paste(myARG, "start", sep = "_")) 

gfa_write(newGfa, "D:/Desktop/testGFA.gfa")

if(sum(str_detect(newGfa$segments$name, "start")) > 1){
  newGfa = newGfa %>% mergeStartSegments(fixLoops = F)
}
  
test = newGfa$segments



test= newGfa$segments %>% filter(str_detect(name, "start")) %>% pull(sequence)
str_extract(test, "^.{1,100}")

str_extract(test, "^.{0,100}")
# gfa_write(newGfa, "D:/Desktop/testGFA.gfa")

test = allPaths[sapply(allPaths, function(x) "23992" %in% x)]
test = allPaths %>% unlist() %>% table()
