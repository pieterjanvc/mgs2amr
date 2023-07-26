#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)

#Databases 
supplementalDB = "D:/Desktop/SupplementalDatabase.db"
resultsDB = "D:/Desktop/results.db"
revisionsDB = "D:/Desktop/revisions.db"

#Limit analysis to the following pipelineIds (or c() for all)
pipelineIds = c()


sDB = dbConnect(dbDriver("SQLite"), supplementalDB)
rDB = dbConnect(dbDriver("SQLite"), resultsDB)
vDB = dbConnect(dbDriver("SQLite"), revisionsDB)

#Get the results
info = tbl(rDB, "pipeline") %>%  collect()

if(length(pipelineIds) == 0){
  pipelineIds = info$pipelineId
} else {
  info = info %>% filter(pipelineId %in% {{pipelineIds}})
}

argInfo = tbl(rDB, "argInfo") %>% 
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
  mutate(grp = overlap(qstart, qend)) %>% ungroup %>% collect()
#----------------

#Get the ground truth
genotypes = tbl(sDB, "genotypes") %>% 
  filter(Run %in% local(info$isolate)) %>% 
  filter(hasLink == 1) %>% 
  # left_join(tbl(myConn, "geneLinks")) %>% 
  collect()
bactInfo = tbl(sDB, "sampleInfo") %>% 
  filter(Run %in% local(info$isolate)) %>% collect() %>% 
  left_join(info %>% select(pipelineId, Run = isolate), by = "Run")
ARG = tbl(sDB, "ARG") %>% collect()

#Get more results...
scaffolds = tbl(rDB, "scaffold") %>% 
  filter(pipelineId %in% local(pipelineIds)) %>% collect()
blastOut = tbl(rDB, "blastOut") %>% 
  filter(pipelineId %in% local(pipelineIds)) %>% collect() %>% 
  mutate(geneId = str_extract(qseqid, "\\d+$") %>% as.integer) %>% 
  left_join(ARG %>% select(geneId, gene), by = "geneId") %>% 
  left_join(bactInfo %>% select(pipelineId, Bacterium), 
            by = c("pipelineId"))

#Get MGS2AMR results
mgs2amrInfo = tbl(vDB, "pipeline") %>% filter(name %in% local(info$name)) %>% 
  filter(statusCode == 5) %>% collect() %>% 
  left_join(info %>% select(revPipeId = pipelineId, name), by = "name")
detectedARG = tbl(vDB, "detectedARG") %>% 
  filter(pipelineId %in% local(mgs2amrInfo$pipelineId)) %>% collect() %>% 
  left_join(ARG %>% select(geneId, gene, subtype), by = "geneId")
annotation = tbl(vDB, "annotation") %>% 
  filter(pipelineId %in% local(mgs2amrInfo$pipelineId)) %>% 
  left_join(tbl(vDB, "bactStrains"), by = "accession") %>% 
  left_join(tbl(vDB, "bactTaxa"), by = "taxid") %>% 
  mutate(isolate = paste(genus, species)) %>% collect() %>% 
  left_join(mgs2amrInfo %>% select(pipelineId, revPipeId), by = "pipelineId") %>% 
  left_join(bactInfo %>% select(revPipeId = pipelineId, Bacterium),
            by = "revPipeId")

dbDisconnect(rDB)
dbDisconnect(sDB)
dbDisconnect(vDB)

#Don't use subtypes for gene matching
genotypes = genotypes %>% select(isolate = Run, subtype = gene) %>% 
  left_join(ARG %>% select(gene, subtype), by = c("subtype")) %>% 
  mutate(gene = ifelse(is.na(gene), subtype, gene)) %>% 
  select(-subtype) %>% distinct()


#------------------ MGS2AMR
test = annotation %>% group_by(pipelineId, geneId, genus, species) %>% 
  filter(fullPath == max(fullPath)) %>% slice(1) %>% 
  group_by(pipelineId, geneId) %>% 
  mutate(
    rank = rank(1/fullPath, ties.method = "min"),
    singleTop = sum(rank == 1) == 1,
    match = isolate == Bacterium
  ) %>% ungroup() %>% filter(match)

test = detectedARG %>% left_join(test %>% 
  select(revPipeId, pipelineId, geneId, extension, fullPath, isolate, rank:match),
    by = c("pipelineId", "geneId")) %>% 
 full_join(genotypes)

#------------- VALIDATION

#Analysis
test1 = argInfo %>% select(pipelineId, gene, coverage, qseqid, gene, bitscore_gene = bitscore, grp) %>% distinct() %>% 
  left_join(info %>% select(pipelineId, isolate), by = "pipelineId") %>% 
  mutate(detection = "XP") %>% 
  group_by(pipelineId, qseqid, grp) %>% 
  #Use both cov and bit score to pick best
  filter(coverage * bitscore_gene == max(coverage * bitscore_gene)) %>% 
  ungroup()

test2 = info %>% select(pipelineId, isolate) %>% 
  left_join(genotypes %>% select(isolate, gene), by = "isolate",
            relationship = "many-to-many") %>% mutate(match = T)

result = test1 %>% full_join(test2, by = c("pipelineId", "isolate", "gene")) %>% 
  mutate(
    detection = replace_na(detection, "FN"),
    detection = case_when(
      detection == "XP" & is.na(match) ~ "XP",
      detection == "XP" & !is.na(match) ~ "TP",
      TRUE ~ "FN"
    )
  ) %>%
  mutate(
    grp = replace_na(grp, 0), 
    coverage = replace_na(coverage, 0),
    qseqid = replace_na(qseqid, "not detected")
  ) %>% 
  group_by(pipelineId, qseqid, grp) %>% filter(coverage == max(coverage)) %>% 
  ungroup()


result = result %>% left_join(
  blastOut %>% 
    select(pipelineId, gene, qseqid, sscinames, Bacterium, rank, 
           bitscore_bact = bitscore, singleTop) %>%
    filter(sscinames == Bacterium) %>% 
    mutate(qseqid = str_remove(qseqid, "_\\d+$")) %>% distinct(), 
  by = c("pipelineId", 'gene', "qseqid"))


result = result %>% filter(!is.na(sscinames) | detection == "FN") %>% 
  #COuld be useful
  # filter(coverage > 0.9 | detection == "FN") %>% 
  arrange(pipelineId, gene, bitscore_gene)

result %>% group_by(pipelineId) %>% summarise(
  nGenes = sum(detection %in% c("TP","FN")), detected = sum(detection == "TP"),
  perc = detected / nGenes
)

### NEED TO FIND A WAY FOR OTHER FP GENES

