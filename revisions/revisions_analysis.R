#/////////////////////////////
# ---- MGS2AMR revisions ----
#////////////////////////////

suppressMessages(library(tidyverse))
library(RSQLite)

#Databases 
supplementalDB = "D:/Desktop/SupplementalDatabase.db"
metaspadesDB = "D:/Desktop/metaspadesPipe.db"
mgs2amrRevDB = "D:/Desktop/mgs2amrRev.db"
monitorDB = "D:/Desktop/monitorDB.db" 

#Limit analysis to the following spPipeIds (or c() for all)
spPipeIds = c(2:6)

sDB = dbConnect(dbDriver("SQLite"), supplementalDB)
spDB = dbConnect(dbDriver("SQLite"), metaspadesDB)
mgsDB = dbConnect(dbDriver("SQLite"), mgs2amrRevDB)

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
# overlap = function(x, y, amount = 50){
#   grp = 1
#   if(length(x) == 1) return(grp)
#   out = rep(1, length(x))
#   z = x
#   x = ifelse(x > y, y, x)
#   y = ifelse(z < y, y, z)
#   x = c(x, x[length(x)])
#   y = c(y, y[length(y)])
#   for(i in 1:(length(x) - 2)){
#     if(min(y[i],y[i+1]) - max(x[i], x[i+1]) >= amount){
#       out[i+1] = grp
#     } else {
#       grp = grp + 1
#       out[i+1] = grp
#     }
#   }
#   return(out)
# }
# argInfo = argInfo %>% group_by(qseqid) %>% 
#   mutate(grp = overlap(qstart, qend)) %>% ungroup %>% collect()
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

mgsResult = mgsResult %>% 
  mutate(x = case_when(
    is.na(runId) ~ "FN",
    is.na(match) ~ "X",
    match & rank == 1 ~ "TP",
    match & rank > 1 ~ "FN",
    # fullPath >= 2500 ~ "TP",
    # match & rank > 1 & fullPath < 2500 ~ "FN",
    !match ~ "FP")) %>% filter(x != "X") %>% group_by(mgsPipeId) %>% 
  summarise(
  nGenes = sum(x %in% c("TP","FN")), detected = sum(x == "TP"),
  perc = detected / nGenes,
  FP = sum(x == "FP")
)


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
  by = c("spPipeId", 'gene', "qseqid")) %>% 
  mutate(detection = ifelse(rank > 1 & detection ==  'TP', "FN", detection))


spResult = spResult %>% filter(!is.na(sscinames) | detection == "FN") %>% 
  #COuld be useful
  # filter(coverage > 0.9 | detection == "FN") %>% 
  arrange(spPipeId, gene, bitscore_gene)

spResult = spResult %>% group_by(spPipeId) %>% summarise(
  nGenes = sum(detection %in% c("TP","FN")), detected = sum(detection == "TP"),
  perc = detected / nGenes
)

### NEED TO FIND A WAY FOR OTHER FP GENES
# utils::browseURL(sprintf("https://www.ncbi.nlm.nih.gov/pathogens/isolates#%s", 
#                          bactInfo$Run[1]))

#------------- MONITORING
mgsDB = dbConnect(dbDriver("SQLite"), mgs2amrRevDB)
monDB = dbConnect(dbDriver("SQLite"), monitorDB)

memZize = function(data, unit = "GB", round = 3){
  
  if(unit == "MB"){
    conversion = c(Mbytes = 1, Gbytes = 1024)
  } else if (unit == "GB"){
    conversion = c(Mbytes = 1/1024, Gbytes = 1)
  } else {
    stop("Unit needs to be 'MB' or 'GB'")
  }
  
  data = lapply(data, str_split, pattern = "\\s+", simplify = T)
  sapply(data, "[[", 1) %>% as.numeric() * 
    unname(conversion[sapply(data, "[[", 2)]) %>% round(round)
}

#Get info from MGS2AMR
info = tbl(mgsDB, "pipeline") %>% select(pipelineId, name) %>% collect()

runIds = tbl(mgsDB, "logs") %>% collect() %>% 
  filter(str_detect(actionName, 
                    paste("pipelineId", info$pipelineId, collapse = "|"))) %>% 
  transmute(runId, pipelineId = str_extract(actionName, "(?<=pipelineId )\\d+") %>% 
              as.integer())
runIds = bind_rows(
  runIds, 
  tbl(mgsDB, "scriptUse") %>% filter(pipelineId %in% local(info$pipelineId)) %>% 
    select(runId, pipelineId) %>% collect()
)

#Only get logs for pipelined that completed in a single run (i.e. no errors)
logs = tbl(mgsDB, "logs") %>% filter(actionId == 1) %>% 
  filter(runId %in% local(runIds$runId)) %>% collect() %>% 
  left_join(runIds, by = "runId") %>% 
  left_join(tbl(mgsDB, "scriptUse") %>% select(runId, start, end) %>% 
              collect(), by = "runId") %>% 
  mutate(across(c(timeStamp, start, end), 
                \(x) as.POSIXct(x, origin = "1970-01-01 UTC"))) %>% 
  group_by(pipelineId) %>% 
  filter(n_distinct(actionName) == 4 & n() == 4) %>% 
  mutate(timeStamp = timeStamp - min(timeStamp)) %>% ungroup()

# logs = tbl(mgsDB, "logs") %>% filter(actionId == 1) %>% left_join(
#   tbl(mgsDB, "scriptUse") %>% filter(status == "finished"), by = "runId"
# ) %>% filter(pipelineId %in% local(c(0,info$pipelineId))) %>% collect() %>% 
#   mutate(
#     pipelineId = ifelse(pipelineId == 0, NA, pipelineId),
#     across(c(timeStamp, start, end), \(x) as.POSIXct(x, origin = "1970-01-01 UTC"))) %>% 
#   # fill(pipelineId) %>% 
#   left_join(x, by = c("runId" = "rId")) %>% mutate(pipelineId = ifelse(is.na(pipelineId), pId, pipelineId)) %>% 
#   group_by(pipelineId) %>% 
#   filter(n_distinct(actionName) == 4 & n() == 4) %>% 
#   mutate(timeStamp = timeStamp - min(timeStamp)) %>% ungroup()

#Get resource monitoring
# x = data.frame(
#   note = c("SRR13338366_ERR2197823_0.03_7", "SRR13302163_ERR2017444_0.028_10",
#     "SRR4065647_ERR3277347_0.098_12", "ERR1557083_ERR2017415_0.041_8",
#     "SRR3170776_ERR3277283_0.05_11"),
#   pipelineId = c(13,14,16,17,18)
# )
monJobs = tbl(monDB, "job") %>% 
  # filter(note %in% local(mgs2amrInfo$name), !is.na(end)) %>% 
  filter(note == "SRR4065647_ERR3277347_0.098_5") %>% 
  collect() %>% 
  # left_join(mgs2amrInfo %>% select(note = name, pipelineId = mgsPipeId), by = "note")
  mutate(pipelineId = c(4,6,7))
monStats = tbl(monDB, "monitor") %>%  
  filter(jobId %in% local(monJobs$jobId)) %>% collect() %>% 
  left_join(monJobs %>% select(jobId, pipelineId), by = "jobId") %>% 
  left_join(logs %>% filter(actionName == "Start local BLASTn") %>% 
              select(pipelineId, start, end), by = "pipelineId") %>% 
  mutate(timeStamp = as.POSIXct(timeStamp),
         mem = memZize(mem),
         mem = ifelse(between(timeStamp - seconds(20), start, end), mem + 150, mem),
         cpuTime = as.numeric(cpuTime),
         cpuTime = replace_na(cpuTime, 0),
         pipelineId = as.factor(pipelineId)) %>% 
  group_by(jobId) %>% mutate(timeStamp = timeStamp - min(timeStamp)) %>% ungroup()
  

# intercept = info %>% filter(name %in% monJobs$note) %>% pull(pipelineId)
intercept = c(4,6,7)
intercept = logs %>% filter(pipelineId %in% intercept) %>% 
  group_by(pipelineId) %>% 
  mutate(actionName = str_remove(actionName, " for pipelineId.*"),
         labelx = timeStamp + (lead(timeStamp) - timeStamp) / 2,
         labelx = ifelse(is.na(labelx), timeStamp + 150, labelx),
         label = LETTERS[1:4]) %>% 
  ungroup()

ggplot(monStats %>% filter(!is.na(pipelineId)) %>% mutate(mem = mem + 0.001), 
       aes(x = timeStamp, y = mem)) + geom_line() +
  geom_vline(data = intercept, aes(xintercept = timeStamp), linetype = "dashed") +
  geom_text(data = intercept, aes(x = labelx, y = 100, label = label), color = "gray50") +
  facet_wrap(vars(pipelineId), scales = "fixed") +
  # scale_y_log10()+
  labs(x = "PIpeline run time (sec)", y = "Memory usage (GB)") +
  theme_minimal()

ggplot(monStats, aes(x = timeStamp, y = cpuTime)) + geom_line() + theme_minimal()

dbDisconnect(mgsDB)
dbDisconnect(monDB)


x = tbl(spDB, "log") %>% filter(pipelineId != 1) %>% collect() %>% 
  arrange(pipelineId, timeStamp) %>% mutate(timeStamp = as.POSIXct(timeStamp)) %>% 
  group_by(pipelineId) %>% 
  summarise(time = difftime(max(timeStamp), min(timeStamp), units = "secs"))

mean(as.integer(x$time))/60

x = tbl(mgsDB, "pipeline") %>% collect() %>% 
 transmute(pipelineId, time = difftime(as.POSIXct(modifiedTimestamp), 
                                       as.POSIXct(startTimestamp), units = "secs")) %>% 
  mutate(h = as.integer(time) / 3600)
mean(x$time)
           