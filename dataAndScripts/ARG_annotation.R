#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse)
library(RSQLite)

tempFolder = "temp/diffMixInPerc0.05_1604442644/"

myConn = dbConnect(SQLite(), "dataAndScripts/meta2amr.db")
ARG = dbReadTable(myConn, "ARG")
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
  list.files(tempFolder, full.names = T,
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
clusterOut = read.table(paste0(tempFolder, "blastSegments.out")) %>% 
  select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
  mutate(clusterId = ifelse(clusterId == "*", segmentId, clusterId)) %>% 
  extract(segmentId, into = c("geneId", "segment"), regex = "^([^_]+)_(.*)", remove = F) %>% 
  filter(clusterId %in% blastOut$query_title) %>% arrange(clusterId)

blastOut = clusterOut %>% 
  left_join(blastOut, by = c("clusterId" = "query_title"))

### Here starts speculation ###

#Use the covermean function so summarise data per gene, bact and segment
blastOut = blastOut %>% group_by(geneId, genus, species, segment, plasmid) %>% 
  summarise(
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

myGene = "2610"
myGFA = gfa_read(sprintf("%sgenesDetected/simplifiedGFA/%s_simplified.gfa", tempFolder, myGene))
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

result = paths[test %>% filter(geneId == myGene) %>% pull(segment)]

result = map_df(result, function(myPath){
  segmentOrder = myPath$segmentOrder[!str_detect(myPath$segmentOrder, "start$")]
  if(length(segmentOrder) > 0){
    data.frame(step = 1:length(segmentOrder), segment = segmentOrder) %>% 
      left_join(test %>% filter(geneId == myGene), by = "segment") %>% 
      left_join(myGFA$segments %>% select(name, LN, KC), by = c("segment" = "name")) %>% 
      mutate(segment = myPath$endSegment) 
  } else {
    data.frame()
  }
})


#------------

