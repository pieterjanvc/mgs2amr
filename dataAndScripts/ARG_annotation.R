#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse)
library(RSQLite)

tempFolder = formatPath("temp/", endWithSlash = T)
sampleName = "testMix2WithBG_1605708377"
minBlastLength = 250

myConn = dbConnect(SQLite(), "dataAndScripts/meta2amr.db")
ARG = dbReadTable(myConn, "ARG")

### FIX mismatch runId in blast prep!! 
genesDetected = dbGetQuery(
  myConn, 
  "SELECT a.runId, a.segmentLength, a.kmerCount, a.n, a.coverage, c.*
  FROM detectedARG a, blastSubmissions as b, ARG as c
  WHERE a.runId == 27 AND b.runId == 28 AND 
  b.tempName == 'testMix2WithBG_1605708377' AND c.geneId = a.geneId") %>% 
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

#Use the covermean function so summarise data per gene, bact and segment
blastOut = blastOut %>% group_by(hitId, geneId, genus, species, segment, plasmid) %>% 
  summarise(
    score_mean = coverMean(query_from, query_to, query_len, score, max),
    bit_score = coverMean(query_from, query_to, query_len, bit_score),
    evalue = coverMean(query_from, query_to, query_len, evalue, 
                       fun = function(x) ifelse(sum(x) == 0, 0, min(x[x != 0]))),
    coverage = coverMean(query_from, query_to, query_len, 1.0),
    n = coverMean(query_from, query_to, query_len, n()),
    query_len = query_len[1], .groups = "drop")

#Get the top for each gene
test = blastOut %>% group_by(geneId, segment) %>% 
  filter(bit_score == max(bit_score)) %>%  filter(coverage > 0.9) 

# test1 = test %>% group_by(geneId, segment) %>% %>% group_by(geneId, genus, species, plasmid) %>% 
#   summarise(bit_score = sum(bit_score), coverage = mean(coverage), n = n(), query_len = sum(query_len)) %>% 
#   arrange(geneId, desc(bit_score)) %>% group_by(geneId) %>%  top_n(2)


#Create a coloured gfa
segmentData = test %>% 
  mutate(info = paste(str_extract(genus, "..."), 
                      str_extract(species, "..."),
                      as.integer(bit_score),
                      collapse = " ")) %>%  
  group_by(geneId, segment) %>% 
  summarise(annot = paste(info, collapse = "; "), .groups = "drop")

myGene = "2096"

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
