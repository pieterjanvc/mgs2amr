#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)

#Databases 
supplementalDB = "D:/Desktop/SupplementalDatabase.db"
resultsDB = "D:/Desktop/results.db"

#Limit analysis to the following pipelineIds (or c() for all)
pipelineIds = c()

#Get the results
myConn = dbConnect(dbDriver("SQLite"), resultsDB)

info = tbl(myConn, "pipeline") %>%  collect()

if(length(pipelineIds) == 0){
  pipelineIds = info$pipelineId
} else {
  info = info %>% filter(pipelineId %in% {{pipelineIds}})
}

argInfo = tbl(myConn, "argInfo") %>% 
  filter(pipelineId %in% local(pipelineIds)) %>% collect()

#---------------- REMOVE NEXT RUN!!!
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
argInfo = argInfo %>% group_by(qseqid) %>% 
  mutate(grp = overlap(qstart, qend)) %>% 
#----------------

scaffolds = tbl(myConn, "scaffold") %>% 
  filter(pipelineId %in% local(pipelineIds)) %>% collect()
blastOut = tbl(myConn, "blastOut") %>% 
  filter(pipelineId %in% local(pipelineIds)) %>% collect() %>% 
  mutate(geneId = str_extract(qseqid, "\\d+$") %>% as.integer) %>% 
  left_join(ARG %>% select(geneId, gene), by = "geneId")

dbDisconnect(myConn)

#Get the ground truth
myConn = dbConnect(dbDriver("SQLite"), supplementalDB)

genotypes = tbl(myConn, "genotypes") %>% 
  filter(Run %in% local(info$isolate)) %>% 
  filter(hasLink == 1) %>% 
  # left_join(tbl(myConn, "geneLinks")) %>% 
  collect()
bactInfo = tbl(myConn, "sampleInfo") %>% 
  filter(Run %in% local(info$isolate)) %>% collect()
ARG = tbl(myConn, "ARG") %>% collect()
dbDisconnect(myConn)

#Don't use subtypes for gene matching
genotypes = genotypes %>% select(isolate = Run, subtype = gene) %>% 
  left_join(ARG %>% select(gene, subtype), by = c("subtype")) %>% 
  mutate(gene = ifelse(is.na(gene), subtype, gene)) %>% 
  select(-subtype) %>% distinct()

#Analysis
test1 = argInfo %>% select(pipelineId, gene, coverage) %>% distinct() %>% 
  left_join(info %>% select(pipelineId, isolate), by = "pipelineId") %>% 
  mutate(detection = "TP")

test2 = info %>% select(pipelineId, isolate) %>% 
  left_join(genotypes %>% select(isolate, gene), by = "isolate",
            relationship = "many-to-many") %>% mutate(match = T)

result = test1 %>% full_join(test2, by = c("pipelineId", "isolate", "gene")) %>% 
  mutate(
    detection = replace_na(detection, "FN"),
    detection = case_when(
      detection == "TP" & is.na(match) ~ "FP",
      detection == "TP" & !is.na(match) ~ "TP",
      TRUE ~ "FN"
    )
  ) %>% select(-match) %>% arrange(pipelineId, gene) %>% 
  left_join(argInfo %>% select(pipelineId, qseqid, gene, grp) %>% distinct(), 
            by = c("pipelineId", "gene")) %>% 
  mutate(
    grp = replace_na(grp, 0), 
    coverage = replace_na(coverage, 0),
    qseqid = replace_na(qseqid, "not detected")
  ) %>% 
  group_by(pipelineId, qseqid, grp) %>% filter(coverage == max(coverage)) %>% 
  ungroup()


result = result %>% left_join(
  blastOut %>% 
    select(pipelineId, gene, qseqid, sscinames, rank, bitscore, singleTop) %>%
    filter(sscinames == bactInfo$Bacterium) %>% 
    mutate(qseqid = str_remove(qseqid, "_\\d+$")) %>% distinct(), 
  by = c("pipelineId", 'gene', "qseqid"))

result = result %>% filter(!is.na(sscinames) | detection == "FN") %>% 
  filter(coverage > 0.9 | detection == "FN")

result %>% group_by(pipelineId) %>% summarise(
  nGenes = sum(detection %in% c("TP","FN")), detected = sum(detection == "TP"),
  perc = detected / nGenes
)


