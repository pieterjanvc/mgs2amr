#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)

#Databases 
supplementalDB = "D:/Desktop/SupplementalDatabase.db"
metaspadesDB = "D:/Desktop/results.db"
mgs2amspDB = "D:/Desktop/revisions.db"

#Limit analysis to the following spPipeIds (or c() for all)
spPipeIds = c()


sDB = dbConnect(dbDriver("SQLite"), supplementalDB)
spDB = dbConnect(dbDriver("SQLite"), metaspadesDB)
mgsDB = dbConnect(dbDriver("SQLite"), mgs2amspDB)

#Get the results
info = tbl(spDB, "pipeline") %>%  collect() %>% 
  rename(spPipeId = pipelineId)

if(length(spPipeIds) == 0){
  spPipeIds = info$spPipeId
} else {
  info = info %>% filter(spPipeId %in% {{spPipeIds}})
}

argInfo = tbl(spDB, "argInfo") %>% 
  filter(pipelineId %in% local(spPipeIds)) %>% collect() %>% 
  rename(spPipeId = pipelineId)

mgs2amrInfo = tbl(mgsDB, "pipeline") %>% filter(name %in% local(info$name)) %>% 
  filter(statusCode == 5) %>% collect() %>% 
  left_join(info %>% select(name, spPipeId), by = "name") %>% 
  rename(mgsPipeId = pipelineId)

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
  collect() %>% rename(isolate = Run) %>% 
  left_join(info %>% select(spPipeId, isolate), by = "isolate",
            relationship = "many-to-many") %>% 
  left_join(mgs2amrInfo %>% select(spPipeId, mgsPipeId), by = "spPipeId")

bactInfo = tbl(sDB, "sampleInfo") %>% 
  filter(Run %in% local(info$isolate)) %>% collect() %>% 
  left_join(info %>% select(spPipeId, Run = isolate), by = "Run") %>% 
  left_join(mgs2amrInfo %>% select(spPipeId, mgsPipeId), by = "spPipeId")

ARG = tbl(sDB, "ARG") %>% collect()

#Get more results...
scaffolds = tbl(spDB, "scaffold") %>% 
  filter(pipelineId %in% local(spPipeIds)) %>% collect() %>% 
  rename(spPipeId = pipelineId)
blastOut = tbl(spDB, "blastOut") %>% 
  filter(pipelineId %in% local(spPipeIds)) %>% collect() %>% 
  mutate(geneId = str_extract(qseqid, "\\d+$") %>% as.integer) %>% 
  left_join(ARG %>% select(geneId, gene), by = "geneId") %>% 
  rename(spPipeId = pipelineId) %>% 
  left_join(bactInfo %>% select(spPipeId, Bacterium),
            by = c("spPipeId"))

#Get MGS2AMR results
# mgs2amrInfo = tbl(mgsDB, "pipeline") %>% filter(name %in% local(info$name)) %>% 
#   filter(statusCode == 5) %>% collect() %>% 
#   left_join(info %>% select(revPipeId = pipelineId, name), by = "name")
detectedARG = tbl(mgsDB, "detectedARG") %>% 
  filter(pipelineId %in% local(mgs2amrInfo$mgsPipeId)) %>% collect() %>% 
  left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>% 
  rename(mgsPipeId = pipelineId)

annotation = tbl(mgsDB, "annotation") %>% 
  filter(pipelineId %in% local(mgs2amrInfo$mgsPipeId)) %>% 
  left_join(tbl(mgsDB, "bactStrains"), by = "accession") %>% 
  left_join(tbl(mgsDB, "bactTaxa"), by = "taxid") %>% 
  mutate(isolate = paste(genus, species)) %>% collect() %>% 
  rename(mgsPipeId = pipelineId) %>% 
  # left_join(mgs2amrInfo %>% select(mgsPipeId, revPipeId), by = "pipelineId") %>% 
  left_join(bactInfo %>% select(mgsPipeId, Bacterium),
            by = "mgsPipeId")

dbDisconnect(spDB)
dbDisconnect(sDB)
dbDisconnect(mgsDB)

#Don't use subtypes for gene matching
genotypes = genotypes %>% 
  select(isolate, subtype = gene, type, core, spPipeId, mgsPipeId) %>% 
  left_join(ARG %>% select(gene, subtype), by = c("subtype")) %>% 
  mutate(gene = ifelse(is.na(gene), subtype, gene)) %>% 
  select(-subtype) %>% distinct()

# test = info %>% select(pipelineId, isolate) %>% 
#   left_join(genotypes %>% select(isolate, gene), by = "isolate",
#             relationship = "many-to-many") %>% mutate(match = T) %>% 
#   left_join(mgs2amrInfo %>% select(pipelineId, revPipeId),
#             by = "pipelineId")


#------------------ MGS2AMR
mgsResult = annotation %>% group_by(mgsPipeId, geneId, genus, species) %>% 
  filter(fullPath == max(fullPath)) %>% slice(1) %>% 
  group_by(mgsPipeId, geneId) %>% 
  mutate(
    rank = rank(1/fullPath, ties.method = "min"),
    singleTop = sum(rank == 1) == 1,
    match = isolate == Bacterium
  ) %>% ungroup() %>% filter(match)

mgsResult = detectedARG %>% left_join(mgsResult %>% 
  select(mgsPipeId, geneId, extension, fullPath, isolate, rank:match),
    by = c("mgsPipeId", "geneId")) %>% 
  full_join(genotypes %>% select(gene, mgsPipeId, core) %>% 
              filter(!is.na(mgsPipeId)), by = c("gene", "mgsPipeId")) %>% 
  mutate(fullPath = replace_na(fullPath, 0)) %>% 
  group_by(mgsPipeId, ARGgroup) %>% filter(fullPath == max(fullPath))

#------------- SPADES VALIDATION

#Analysis
test1 = argInfo %>% select(spPipeId, gene, coverage, qseqid, gene, 
                           bitscore_gene = bitscore, grp) %>% distinct() %>% 
  left_join(info %>% select(spPipeId, isolate), by = "spPipeId") %>% 
  mutate(detection = "XP") %>% 
  group_by(spPipeId, qseqid, grp) %>% 
  #Use both cov and bit score to pick best
  filter(coverage * bitscore_gene == max(coverage * bitscore_gene)) %>% 
  ungroup()

# test2 = info %>% select(pipelineId, isolate) %>% 
#   left_join(genotypes %>% select(isolate, gene), by = "isolate",
#             relationship = "many-to-many") %>% mutate(match = T)

spResult = test1 %>% full_join(genotypes %>% mutate(match = T), 
                             by = c("spPipeId", "isolate", "gene")) %>% 
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
  group_by(spPipeId, qseqid, grp) %>% filter(coverage == max(coverage)) %>% 
  ungroup()


spResult = spResult %>% left_join(
  blastOut %>% 
    select(spPipeId, gene, qseqid, sscinames, Bacterium, rank, 
           bitscore_bact = bitscore, singleTop) %>%
    filter(sscinames == Bacterium) %>% 
    mutate(qseqid = str_remove(qseqid, "_\\d+$")) %>% distinct(), 
  by = c("spPipeId", 'gene', "qseqid"))


spResult = spResult %>% filter(!is.na(sscinames) | detection == "FN") %>% 
  #COuld be useful
  # filter(coverage > 0.9 | detection == "FN") %>% 
  arrange(spPipeId, gene, bitscore_gene)

spResult %>% group_by(spPipeId) %>% summarise(
  nGenes = sum(detection %in% c("TP","FN")), detected = sum(detection == "TP"),
  perc = detected / nGenes
)

### NEED TO FIND A WAY FOR OTHER FP GENES
utils::browseURL(sprintf("https://www.ncbi.nlm.nih.gov/pathogens/isolates#%s", 
                         bactInfo$Run[1]))
