#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse)
library(RSQLite)

tempFolder = formatPath("temp/", endWithSlash = T)
sampleName = "testMix2cWithBG_1605805981"
minBlastLength = 250

myConn = dbConnect(SQLite(), "dataAndScripts/meta2amr.db")
ARG = dbReadTable(myConn, "ARG")

### FIX mismatch runId in blast prep!! 
genesDetected = dbGetQuery(
  myConn, 
  "SELECT a.runId, a.segmentLength, a.kmerCount, a.n, a.coverage, c.*
  FROM detectedARG a, ARG as c
  WHERE a.runId == 40 AND c.geneId = a.geneId") %>% 
  distinct()
dbDisconnect(myConn)


# ---- FUNCTIONS ----
#********************
#Get the mean based on the max (default) at each position in the coverage
coverMean = function(from, to, n, val, fun = max){
  apply(
    mapply(function(from, to, n, val, fun){
      x = rep(0.0, n)
      x[from:to] = val
      return(x)
    }, from, to, n, val),
    1, fun
  ) %>% mean()
}


# ---- FILTERING DATA ----
#*************************

#Read all blast output
blastOut = lapply(
  list.files(paste0(tempFolder, sampleName), full.names = T,
             pattern = "blastSegmentsClustered\\d+.json.gz"),
  blast_readResults, outFormat = "dataFrame1") %>% bind_rows()

#Extract data we need 
blastOut = blastOut %>% 
  select(query_title, hitId, taxid, bact = title, bit_score, 
         score, evalue, identity, query_len, query_from, 
         query_to, hit_from, hit_to, align_len) %>% 
  mutate(bact = str_remove_all(bact, "[^\\w\\s]")) %>% 
  extract(bact, c("genus", "species", "extra"), 
          regex = "(\\w+)\\s+(\\w+)\\s+(.*)") %>% 
  filter(!species %in% c("bacterium", "sp")) %>% 
  mutate(plasmid = str_detect(extra, "plasmid|Plasmid"))

#Expand results from clustering segments before blast
clusterOut = read.table(paste0(tempFolder, sampleName, "/blastSegments.out")) %>% 
  select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
  mutate(clusterId = ifelse(clusterId == "*", segmentId, clusterId)) %>% 
  extract(segmentId, into = c("geneId", "segment"), regex = "^([^_]+)_(.*)", remove = F) %>% 
  filter(clusterId %in% blastOut$query_title) %>% arrange(clusterId)

blastOut = clusterOut %>% 
  left_join(blastOut, by = c("clusterId" = "query_title"))


#It seems that using multiple matches within the same genome does not make sense,
#since it would not have been split if it would be extended with a gap.
#So better get rid of all but the best match? 
# We also have to work by hit Id before we start merging
#  group_by(clusterId, hitId)

### Here starts speculation ###
#group_by(hitId, geneId, genus, species, segment, plasmid)
test = blastOut %>% 
  mutate(coverage = align_len / query_len) %>% 
  select(-c(identity:align_len)) %>% 
  group_by(clusterId, hitId) %>% 
  filter(score == max(score)) 

#Create a coloured gfa
segmentData = test %>% 
  mutate(info = paste(str_extract(genus, "..."), 
                      str_extract(species, "..."),
                      as.integer(bit_score),
                      collapse = " ")) %>%  
  group_by(geneId, segment) %>% 
  summarise(annot = paste(info, collapse = "; "), .groups = "drop")

myGene = "5123"

myGFA = gfa_read(sprintf("%s%s/genesDetected/simplifiedGFA/%s_simplified.gfa", 
                         tempFolder, sampleName, myGene))
myGFA = gfa_annotation(
  myGFA, 
  segments = segmentData %>% filter(geneId == myGene) %>% pull(segment), 
  label = segmentData %>% filter(geneId == myGene) %>% pull(annot))

# gfa_write(myGFA, "D:/Desktop/testGFA.gfa")

#Use the path table to see if you can reconstruct path for bacterium all the way to start
start = myGFA$segments %>% filter(str_detect(name, "_start$")) %>% 
  filter(LN == max(LN)) %>% pull(name)
paths = gfa_pathsToSegment(myGFA, start, returnList = T)
names(paths) = sapply(paths, "[", "endSegment")

# result = paths[test %>% filter(geneId == myGene) %>% pull(segment)]

result = map_df(paths, function(myPath){
  # segmentOrder = myPath$segmentOrder[!str_detect(myPath$segmentOrder, "start$")]
  segmentOrder = myPath$segmentOrder
  if(length(segmentOrder) > 0){
    data.frame(step = 1:length(segmentOrder), segment = segmentOrder) %>% 
      left_join(test %>% 
                  filter(geneId == myGene), by = "segment") %>% 
      left_join(myGFA$segments %>% 
                  select(name, LN, KC), by = c("segment" = "name")) %>% 
      mutate(endSegment = myPath$endSegment,
             pathLength = myPath$pathLength)
  } else {
    data.frame()
  }
}) 


test1 = result %>% filter(!is.na(geneId)) %>% 
  group_by(endSegment, taxid, genus, species, extra) %>% 
  mutate(n = n()) %>% 
  group_by(endSegment) %>% filter(n == max(n)) %>% 
  group_by(taxid, genus, species, extra) %>%  
  summarise(score = sum(score), bit_score = sum(bit_score)) %>% 
  group_by(taxid) %>% filter(score == max(score)) %>% 
  arrange(desc(score))


# %>% group_by(endSegment) %>% 
#   mutate(pathLength = pathLength[1] - sum(LN[is.na(geneId)]) + 60*n() - 60)

#Filter by minBlastLength
# test1 = result %>% 
#   # mutate(LN = LN - 30) %>% 
#   # filter(LN >= minBlastLength) %>%
#   group_by(endSegment, genus, species) %>% 
#   mutate(pathLength = pathLength - sum(LN[is.na(geneId)])) %>% 
#   summarise(bit_score = sum(bit_score), n = n(), 
#             pathAssigned = sum(LN) / pathLength[1],
#             .groups = "drop")


test1 = result %>% filter(!is.na(geneId)) %>% 
  mutate(unique = ifelse(segment == endSegment, LN / pathLength, 0)) %>% 
  group_by(endSegment, genus, species, plasmid) %>% mutate(n = n()) %>% 
  group_by(endSegment, plasmid) %>% filter(n == max(n)) %>% 
  group_by(endSegment, genus, species, plasmid)  %>% 
  summarise(bit_score = sum(bit_score), 
            unique = max(unique), 
            quality = weighted.mean(coverage, LN) * pathLength,
            .groups = "drop") %>% 
  filter(unique > 0)

genesDetected %>% filter(geneId == myGene)
