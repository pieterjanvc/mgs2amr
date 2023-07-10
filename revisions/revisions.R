#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

library(tidyverse)
library(RSQLite)

isolate="SRR3345827"
background="ERR3277283"

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

diaData = read.csv("D:/Desktop/diamondOutput.csv", sep = "\t", header = F) 
colnames(diaData) = c("qseqid", "sseqid", "qlen", "slen", "pident", "nident",
                      "length", "evalue", "bitscore", "qstart", "qend")

myConn = dbConnect(dbDriver("SQLite"), "D:/GoogleDrive/My Drive/Varia/UC/SupplementalDatabase.db")
ARG = tbl(myConn, "ARG") %>% collect()
bactInfo = tbl(myConn, "genotypes") %>% filter(Run == isolate) %>% 
  left_join(tbl(myConn, "geneLinks"), by = "gene") %>% collect() %>% 
  filter(!is.na(geneId))
dbDisconnect(myConn)

argInfo = diaData %>% 
  left_join(ARG %>% select(geneId, prot, gene, subtype), by = c("sseqid" = "prot")) %>% 
  mutate(coverage = nident / slen) %>% filter(coverage >= 0.5) %>%
  left_join(bactInfo %>% select(geneId, type), by = "geneId") %>% 
  filter(!is.na(gene)) %>% 
  group_by(qseqid, qstart, qend) %>% filter(evalue == min(evalue)) %>% 
  group_by(qseqid) %>% arrange(qstart, qend) %>% 
  # mutate(overlap = ifelse(qend < lag(qend, default = 1), qend, lag(qend, default = 0)) - 
  #          qstart)
  mutate(grp = overlap(qstart, qend)) %>% 
  group_by(qseqid, grp) %>% filter(bitscore == max(bitscore)) %>% ungroup()

#Link the qsedId (scaffold) to the nodes in the assembly graph
allPaths = readLines("D:/Desktop/scaffolds.paths")
idx = str_detect(allPaths, "^NODE") %>% which()
nodeNames = allPaths[idx]
x = nodeNames %in% argInfo$qseqid %>% which()

allPaths = mapply(function(start, end){
  nodes = allPaths[(start + 1):(end-1)] %>% paste(collapse = "") %>% 
    str_extract_all("\\d+")
  nodes[[1]]
}, start = idx[x], end = idx[x+1])

names(allPaths) = nodeNames[x]
