# detach("package:gfaTools", unload=TRUE)
# library(gfaTools)

#Get all ids (temp)
processed = readLines("/mnt/meta2amrData/pipelineTest/after1200/rerunBlastn.out") %>% 
  str_extract("(?<=pipelineId\\s)\\d+")
processed = processed[!is.na(processed)]
processed = c(2:10, processed)

baseFolder = "/mnt/meta2amrData/meta2amr/"
database = "/mnt/meta2amrData/pipelineTest/after1200/meta2amr.db"
# database = "/mnt/meta2amrData/pipelineTest/extraTests/extraTests.db"
runId = 0
verbose = 1
pipelineIds = NULL
generateReport = F
forceRedo = F
minBlastLength = 250
outfmt = "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart send bitscore score length pident nident qcovs qcovhsp"


myId = "86" #452, 845, 18

# registerDoParallel(cores=5)
# test = foreach(myId = processed[251:400]) %dopar% {
#   
# tryCatch({
  
# test = map_df(as.character(2:10), function(myId){


#FUNCTIONS
softmax = function(vals, normalise = F, log = T){
  if(normalise) vals = vals / max(vals)
  if(log) vals = log(vals)
  return(exp(vals) / sum(exp(vals)))
}

blast_readOutput = function(file, outfmt, separate = T, includeIssues = F, verbose = 1){
  
  #Load the output csv (zipped)
  blastOut = read.table(file, sep = "\t", quote = "", comment.char = "")
  colnames(blastOut) = strsplit(outfmt, " ")[[1]][-1]
  
  #split multiple matches
  blastOut = blastOut %>% 
    mutate(x = str_count(staxids, ";"), y = str_count(sallacc, ";"), 
           z = x == y | x == 0, plasmid = str_detect(salltitles, "plasmid|Plasmid")) 
  
  if(!includeIssues){
    blastOut = blastOut %>% filter(!(x > 0 & !z))
    issues = data.frame()
  } else {
    #Multiple taxid but not the same number of accessions
    issues = blastOut %>% filter(x > 0 & !z) 
  }
  
  if(nrow(issues) > 0){
    warning(nrow(issues), " rows contain ambiguous accession / taxid results")
  }
  
  if(separate){
    return(bind_rows(
      
      #No merged data
      blastOut %>% filter(y == 0),
      #Identical number of accession and taxids 
      blastOut %>% 
        filter(x > 0 & z) %>% 
        separate_rows(sallacc, staxids, sscinames, sep = ";"),
      #One taxId for multiple accessions
      blastOut %>% 
        filter(x == 0 & y > 0 & z) %>% 
        separate_rows(sallacc, staxids, sscinames, sep = ";"),
      issues
      
    ) %>% select(-x, -y, -z) %>% mutate(staxids = as.character(staxids)))
  } else {
    return(blastOut %>% select(-x, -y, -z) %>% 
             mutate(staxids = as.character(staxids))) 
  }
  
}


# test = map_df(as.character(10), function(myId){
myConn = dbConnect(SQLite(), database)

ARG = dbReadTable(myConn, "ARG") 

#Generate a list out of the settings file
settings = readLines(paste0(baseFolder, "settings.txt"))
settings = settings[str_detect(settings,"^\\s*[^=#]+=.*$")]
settings = str_match(settings, "\\s*([^=\\s]+)\\s*=\\s*(.*)")
settings = setNames(str_trim(settings[,3]), settings[,2])

#Get the AMR prediction models
bactGenomeSize = 
  read.csv(sprintf("%sdataAndScripts/bactGenomeSize.csv",baseFolder))

toProcess = dbReadTable(myConn, "pipeline") %>% 
  filter(statusCode > 2)

sampleIndex = which(toProcess$pipelineId == myId)
dbDisconnect(myConn)

sample = toProcess$tempFolder[sampleIndex]
myPipelineId = toProcess$pipelineId[sampleIndex]
sampleName = str_extract(sample, "[^\\/]+(?=_\\d+$)")

inputfileBP = grep("sequences length",
                   readLines(sprintf("%s/metacherchant_logs/log", sample)),
                   value = T, fixed = T) %>%
  str_extract("([\\d'])+(?=\\s\\()") %>%
  str_replace_all("'", "") %>% as.numeric()

#Grab the detected ARG from the previous step
myConn = dbConnect(SQLite(), database,synchronous = NULL)
sqliteSetBusyHandler(myConn, 30000)

genesDetected = tbl(myConn, "detectedARG") %>% 
  filter(pipelineId == myPipelineId) %>% as.data.frame() %>% 
  mutate(cover = ifelse(type == 'noFragments', cover1, cover2))
dbDisconnect(myConn)

newLogs = data.frame(
  timeStamp = as.integer(Sys.time()), 
  actionId = 1, actionName = "Start Annotation & Prediction")

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    sprintf("Start Annotation & Prediction for pipelineId %i (%i/%i)\n", 
                            myPipelineId, sampleIndex, nrow(toProcess)))}

sample

# ---- LOAD DATA ----
#*************************
if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Loading BLASTn output ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                              "Start loading BLASTn output"))

#Read all blast output
blastOut = lapply(
  list.files(sample, full.names = T,
             pattern = "blastSegmentsClustered\\d+.csv.gz"),
  blast_readOutput, outfmt = outfmt) %>% bind_rows() %>% 
  mutate(geneId = str_extract(qseqid, "^([^_]+)"))


expandBlast = list.files(sample, full.names = T, pattern = "expand_\\d+.csv.gz")

if(length(expandBlast) > 0){
  
  #Get the new results
  blastOut2 = lapply(
    list.files(sample, full.names = T,
               pattern = "expand_\\d+.csv.gz"),
    blast_readOutput, outfmt = outfmt) %>% 
    bind_rows() %>% 
    mutate(geneId = str_extract(qseqid, "^([^_]+)"))
  
  blastOut = bind_rows(
    blastOut %>% filter(!qseqid %in% unique(blastOut2$qseqid)),
    blastOut2
  )
  
  rm(blastOut2)
}

#FIX FROM HERE
# backup1 = blastOut 
# blastOut = backup1

#Extract data we need + transform
blastOut = blastOut %>% 
  select(query_title = qseqid, taxid = staxids, accession = sallacc, 
         bact = sscinames, plasmid, bit_score = bitscore, score, 
         identity = nident, query_len = qlen, query_from = qstart, query_to = qend, 
         hit_from = sstart, hit_to = send, align_len = length, geneId) %>% 
  mutate(bact = str_remove_all(bact, "[^\\w\\s]"),
         geneId = as.integer(geneId), taxid = as.integer(taxid)) %>% 
  extract(bact, c("genus", "species", "extra"), 
          regex = "(\\w+)\\s+(\\w+)($|\\s+.*)") %>% 
  #Sp. will be pasted with taxid to make it unique
  mutate(species = ifelse(species == "sp", paste0(species, taxid), species)) %>% 
  filter(!genus %in% c("uncultured", "Uncultured", "mixed", "Bacterium"),
         !species %in% c("bacterium", "Bacterium") & 
           !is.na(species)) %>% 
  #Extract subspecies, strain and plasmid info
  mutate(
    subspecies = str_extract(extra, "(?<=subsp\\s)[^\\s]+"),
    strain = str_extract(extra, "(?<=strain\\s)[^\\s]+") %>% 
      str_remove(" chromosome.*| plasmid.*| complete.*")
  )

#Expand results from clustering segments before blast
clusterOut = read.table(paste0(sample, "/blastSegments.out")) %>% 
  select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
  mutate(
    clusterId = ifelse(clusterId == "*", segmentId, clusterId)
  ) %>% 
  filter(clusterId %in% blastOut$query_title) 

blastOut = clusterOut %>% 
  left_join(blastOut, by = c("clusterId" = "query_title")) %>% 
  #Get the coverage and identity of each alignment
  rowwise() %>% 
  mutate(
    ident = identity / align_len,
    coverage = min(align_len / query_len, 1),
  ) %>% ungroup() 

#Fix some known issues with naming
blastOut[blastOut$genus == "Enterobacter" & 
           blastOut$species == "aerogenes","genus"] = "Klebsiella"


# test = blastOut %>% filter(!segmentId %in% blastSegments$name)
# test = blastSegments %>% filter(!name %in% blastOut$segmentId)

#--
# backup2 = blastOut
# blastOut = backup2
#--

if(verbose > 0){cat("done\n")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                              "Finished loading BLASTn output"))

#Check
blastSegments = 
  read.csv(paste0(sample, "/blastSegments.csv")) %>% 
  select(name, LN, KC) %>% 
  mutate(
    start = ifelse(str_detect(name, "_start$"), T, F),
    geneId = str_extract(name, "\\d+")) %>% 
  #Correct KC wher depth is an outlier
  group_by(geneId) %>% mutate(
    depth = KC / LN,
    z = (depth - mean(depth)) / sd(depth),
    depth = ifelse(z > 1.96,  sum(KC[z <= 1.96 & LN >= 250]) / 
                     sum(LN[z <= 1.96 & LN >= 250]), depth),
    KC = depth * LN
  ) %>% ungroup() %>% select(-geneId, -z, -depth)

blastOut = blastOut %>% 
  left_join(blastSegments, by = c("segmentId" = "name")) %>% 
  mutate(KC = KC * coverage * ident)

# test = blastSegments %>% mutate(
#   found = blastSegments$name %in% blastOut$segmentId) %>% 
#   group_by(found) %>% 
#   summarise(n = n()) %>% mutate(pipelineId = myId, sample = sample)

# return(test)
# })


# if(test[test$found == T, "n"] / sum(test$n) < 0.90) stop("BAD BLAST")


# ---- MAPPING DATA TO GFA ----
#******************************

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Mapping alignments to GFA ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                              "Start mapping"))


segmentsOfInterest = read.csv(paste0(sample, "/segmentsOfInterest.csv"))

#Get information about the blasted seq where they are in the path of GFA
pathData = list.files(
  paste0(sample, "/genesDetected/simplifiedGFA"), 
  pattern = ".gfa", full.names = T)

pathData = map_df(pathData, function(myGFA){
  
  #Read the GFA
  fullGFA = gfa_read(myGFA)
  
  myGene = str_extract(fullGFA$segments$name[1], "^\\d+")
  
  segmentOfInterest = segmentsOfInterest %>% 
    filter(geneId == myGene, type != "fragmentsOnly") 
  
  maxPathId = 0
  result = data.frame()
  
  for(i in 1:nrow(segmentOfInterest)){
    
    myGroup = segmentOfInterest %>% filter(name == segmentOfInterest$name[i]) %>% 
      pull(group)
    
    gfa = gfa_filterSegments(fullGFA, segments = (fullGFA$segments %>% 
                               filter(group == myGroup) %>% pull(name)))
    
    
    #If there are linked segments ...
    if(nrow(gfa$links) > 0){
      
      #Get all semenents in path to start
      x = gfa_pathsToSegment(gfa, segmentOfInterest$name[i], returnList = T, 
                             pathSegmentsOnly = T, verbose = F, allowReverse = F,
                             maxDistance = 2000) %>% 
        map_df(function(path){
          data.frame(
            pathId = path$id,
            orientation = path$orientation,
            pathType = path$pathType,
            segmentId = path$segmentOrder,
            dist = path$dist,
            LN = path$LN,
            KC = path$KC)
        }) %>% 
        mutate(geneId = str_extract(segmentId, "^\\d+"),
               group = myGroup)
      
      
      #If no paths, return empty frame
      if(nrow(x) != 0 ){
        
        x = x %>% 
          group_by(geneId, pathId) %>% mutate(
            pathId = pathId + maxPathId,
            order = n():1,
            #These files contain full GFA structures, no fragments
            type = segmentOfInterest$type[i]) %>% 
          ungroup()
        
        maxPathId = max(x$pathId)
        
        result = bind_rows(result, x)
        
      }
      
      
      
    } else if(nrow(gfa$segments) > 0){
      
      #Case where there are no links, just one or more segments
      x = gfa$segments %>% 
        select(segmentId = name, KC, LN) %>% 
        mutate(
          pathId = 1:n() + maxPathId,
          orientation = 0,
          geneId = str_extract(segmentId, "^\\d+"),
          dist = ifelse(str_detect(segmentId, "_start$"), -LN, 0),
          order = 1,
          type = "full"
        )
      
      maxPathId = max(x$pathId)
      
      result = bind_rows(result, x)
      
    } 
  }
  
  return(result)
  
}) %>% mutate(geneId = as.integer(geneId))


#Calculate the score to bit_score conversion factor for each gene
extraBits = blastOut %>% filter(ident == 1, coverage == 1) %>% 
  group_by(geneId) %>% 
  summarise(bitConst = mean(bit_score / score), .groups = "drop")

extraBits = pathData %>% 
  filter((LN < minBlastLength & dist >= 0 ) | (dist < 0 & LN < 75)) %>% 
  group_by(geneId, pathId) %>% 
  summarise(score = sum(LN) - 30*n(), .groups = "drop") %>% 
  left_join(extraBits, by = "geneId") %>% 
  mutate(bitScore = bitConst * score)

#Only work with segments >= 250 in paths
pathData = pathData %>% 
  filter(LN >= minBlastLength | (dist < 0 & LN > 74)) %>% 
  group_by(geneId, pathId) %>% 
  arrange(geneId, pathId, desc(order)) %>% 
  mutate(order = n():1) %>% 
  ungroup()

#Also add the fragment data to the paths
fragments = gfa_read(paste0(sample, "/fragmentGFA.gfa"))

if(nrow(fragments$segments) > 0){
  #Fragments were merged, so treated as a non-start piece
  fragments = fragments$segments %>% 
    select(segmentId = name, LN, KC, geneId) %>% 
    group_by(geneId) %>% 
    mutate(pathId = 1, orientation = 0:(n()-1),
           order = 1, type = "fragment", dist = -1,
           geneId = as.integer(geneId), depth = KC / LN)
} else {
  fragments = data.frame()
}

pathData = bind_rows(pathData, fragments)  %>% group_by(geneId) %>% 
  mutate(orientation = ifelse(type == "fragment" & any(type == "full"), 
                              abs(orientation[type == "full"][1]-1), orientation)) %>% 
  ungroup()


pathData = pathData %>% group_by(geneId, pathId, orientation) %>% 
  mutate(totalLN = sum(LN[dist >= 0]) - 30*sum(dist >= 1)) %>% ungroup()

if(verbose > 0){cat("done\n")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, 
                              "Finished mapping"))


#Add the grouping info to detected genes
genesDetected = genesDetected %>%
  left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>%
  select(pipelineId, geneId, gene, subtype, everything()) %>% 
  mutate(ARGgroup = 1:n())

# ---- Filter / group bacteria ----
#**********************************

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Group overlapping bacterial species ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, 
                              "Start grouping bacteria"))


#Make sure that taxid always has the same genus / species name and vice versa
blastOut = blastOut %>% 
  mutate(bact = paste(genus, species)) %>%
  group_by(taxid, genus, species) %>% 
  mutate(myCount = n()) %>% 
  group_by(taxid) %>% 
  mutate(genus = genus[myCount == max(myCount)][1],
         species = species[myCount == max(myCount)][1]) %>% 
  group_by(taxid, genus, species) %>% 
  mutate(myCount = n()) %>% 
  group_by(genus, species) %>% 
  mutate(taxid = taxid[myCount == max(myCount)][1]) %>% 
  ungroup()

blastOut$start[blastOut$start & ! blastOut$segmentId %in% 
                 segmentsOfInterest$name[segmentsOfInterest$type != "fragmentsOnly"]] = F

#Add the group (membership of subgraphs)
blastOut = blastOut %>% left_join(
  pathData %>% select(segmentId, group) %>% distinct(), by = "segmentId") %>% 
  mutate(geneId = as.integer(geneId))

# When a strain flanks both sides of a the start segment, but the segment itself has no match,
# this likely results from too strict cutoffs by BLASTn. Add a start match with 99%
# of the bit_score of the other matches to keep the accession in the running
addSeedHit = pathData %>% 
  select(geneId, orientation, segmentId, order, dist, group) %>% 
  group_by(geneId, orientation) %>% 
  filter(order < 3) %>% distinct() %>% ungroup() %>% 
  left_join(
    blastOut %>% 
      select(segmentId, accession, genus, species, hit_from, hit_to, 
             start, bit_score, identity, coverage, LN, KC) %>% 
      distinct(),
    by = "segmentId"
  ) %>% group_by(geneId, accession) %>% 
  filter(n_distinct(orientation[order > 1]) > 1) 

#Get the bit_scores of the start segment
otherSeedHits = addSeedHit %>% ungroup() %>% filter(any(order == 1)) %>% 
  select(-c(dist, group, orientation, accession:hit_to)) %>% 
  mutate(geneId = as.integer(geneId)) %>% 
  filter(order == 1) %>% group_by(geneId) %>% 
  filter(bit_score == max(bit_score)) %>% distinct() %>% ungroup() %>% 
  select(-order) %>% mutate(bit_score = 0.99 * bit_score)

# Create the start hit in the genome by looking at the positions of the matches on 
# either side of the start segment, and 
addSeedHit = addSeedHit %>% select(geneId:accession, hit_from, hit_to) %>% 
  filter(all(order != 1)) %>%
  mutate(x = hit_from == max(c(hit_from, hit_to)) | 
           hit_to == max(c(hit_from, hit_to))) %>% 
  group_by(geneId, accession, orientation) %>% 
  mutate(x = ifelse(any(x == T), T, F)) %>% 
  ungroup() %>% rowwise() %>% 
  mutate(
    hit_from = ifelse(x, hit_from - dist, hit_from + dist),
    hit_to = ifelse(x, hit_to - dist, hit_to + dist),
    d1 = ifelse(x, min(hit_from, hit_to), max(hit_from, hit_to)),
         geneId = as.integer(geneId)) %>% 
  group_by(geneId, accession) %>%
  mutate(d2 = mean(d1[x]) - mean(d1[!x])) %>% ungroup() %>% 
  left_join(ARG %>% select(geneId, nBases), by = "geneId") %>% 
  mutate(d3 = abs(d2 - nBases) - 30) %>% 
  filter(!is.na(d3 >= 0)) %>% filter(d3 <= 250) %>% 
  group_by(geneId, accession, group, nBases) %>% 
  summarise(    d1 = min(d1), d2 = max(d1), .groups = "drop") %>% mutate(
    hit_from = d1 + ((d1 - d2) - nBases) / 2,
    hit_to = d2 - ((d1 - d2) - nBases) / 2
  )

addSeedHit = otherSeedHits %>% filter(geneId %in% addSeedHit$geneId) %>% 
  left_join(addSeedHit %>% select(geneId, accession, hit_from, hit_to, group)) %>% 
  left_join(
    blastOut %>% 
      select(accession, taxid, genus, species, extra, 
             subspecies, plasmid) %>% 
      distinct(), by  = "accession")

#Add the new start 'matches' to the blast output
blastOut = bind_rows(blastOut, addSeedHit)

#Only consider segments in the graph that are also the same distance
#in the actual genome alignments

#  IDEA : 3083_unitig229_start. If accession is present on both sides,
# it automatically passes the test, even if not in start (slightly lower match)
#
tempVar = blastOut %>% 
  # filter(geneId == "5166") %>%
  group_by(geneId, accession, group) %>% 
  filter(any(start)) %>% 
  mutate(sMin = ifelse(hit_from < hit_to, hit_from, hit_to),
         sMax = ifelse(hit_from >= hit_to, hit_from, hit_to),
         after = ifelse(sMin > sMin[start], T, F),
         dist = ifelse(after, sMin - sMax[start], sMin[start] - sMax) + 29) %>%
  ungroup() %>% 
  left_join(pathData %>% select(segmentId, dist1 = dist) %>% 
              distinct(), by = "segmentId") %>% 
  rowwise() %>% 
  #Repeats can make segments shorter because loop is cut off, so provide some slack
  mutate(correct = between(dist, dist1 -500, dist1 + 500) | start) %>% 
  filter(correct)

wrongLocation = blastOut %>% filter(!segmentId %in% tempVar$segmentId)
# tempVar = blastOut %>% filter(!segmentId %in% wrongLocation$segmentId)
blastOut = bind_rows(
  tempVar,
  blastOut %>% 
    filter(geneId %in% (genesDetected %>% 
                          filter(!geneId %in% segmentsOfInterest$geneId) %>% 
                          pull(geneId)))
  )
rm(tempVar)

allBact = blastOut %>%
  # filter(geneId %in% c("4793")) %>%
  select(segmentId, geneId, bit_score, coverage,
         accession, taxid, genus, species, extra, plasmid, KC, LN, start) %>%
  group_by(segmentId, geneId, accession) %>%
  filter(bit_score == max(bit_score)) %>% 
  dplyr::slice(1) %>% ungroup() %>%
  mutate(pathScore = bit_score * coverage) %>%
  # mutate(pathScore = bit_score) %>% 
  group_by(segmentId, accession, plasmid) %>% 
  filter(pathScore == max(pathScore)) %>% dplyr::slice(1) %>% ungroup() %>% 
  left_join(pathData %>% 
              select(pathId, segmentId, order, orientation) %>% 
              mutate(orientation = ifelse(str_detect(segmentId, "_start$"),
                                          -1, orientation)) %>% 
              distinct(), by = "segmentId") %>% 
  filter(!is.na(order))  %>% 
  arrange(geneId, pathId, accession, plasmid, order) %>% 
  group_by(geneId, pathId, accession, plasmid) %>% 
  filter(order == 1:n()) %>% 
  
  
  # group_by(geneId, pathId, accession, plasmid, orientation) %>% 
  # filter(all(2:max(order) %in% order) | all(order == 1)) %>% 
  
  
  # group_by(segmentId, accession) %>% filter(pathScore == max(pathScore)) %>% 
  # dplyr::slice(1) %>% group_by(geneId, accession, pathId) %>% 
  # mutate(fullPath = sum(pathScore))  
  
  #Calculate path score
  group_by(geneId, accession, plasmid, orientation, pathId) %>%
  mutate(fullPath = sum(pathScore)) %>%
  #Pick best path for each orientation
  group_by(geneId, accession, plasmid, orientation) %>%
  filter(pathId == pathId[fullPath == max(fullPath)][1]) %>%
  #Sum the 3 orientations to get full path score
  group_by(geneId, accession, plasmid) %>%
  mutate(fullPath = sum(pathScore)) %>%
  group_by(geneId, accession) %>%
  filter(plasmid == plasmid[fullPath == max(fullPath)][1]) %>% 
  ungroup()

#Calculate the exta bits to add to the total path score
test = allBact %>% filter(orientation != -1) %>% 
  select(geneId, pathId, accession) %>% distinct() %>% 
  left_join(extraBits, by = c("geneId", "pathId")) %>%
  group_by(geneId, accession) %>% 
  mutate(score = sum(score, na.rm = T),
         bitScore = sum(bitScore, na.rm = T)) %>% 
  group_by(geneId, accession) %>% 
  summarise(score = max(score), bitConst = max(0, bitConst, na.rm = T),
                        bitScore = max(bitScore), .groups = "drop")
#Add the extra bits
allBact = allBact %>% 
  left_join(test, by = c("geneId", "accession")) %>%
  mutate(score = replace_na(score, 0), bitScore = replace_na(bitScore, 0)) %>% 
  mutate(fullPath = fullPath + bitScore)
  
allBact = allBact %>% 
  left_join(
    allBact %>% group_by(geneId, accession) %>% 
      summarise(fullPath = max(fullPath), .groups = "drop") %>% 
      group_by(accession) %>% 
      summarise(fullPath = sum(fullPath),.groups = "drop") %>% 
      transmute(accession, rank = rank(fullPath)),
    by = "accession"
  ) %>% 
  group_by(geneId, accession) %>% 
  mutate(path0 = pathId[orientation == 0][1],
         maxOrder0 = max(0,order[orientation == 0]),
         path1 = pathId[orientation == 1][1],
         maxOrder1 = max(0,order[orientation == 1]),
         LN = sum(LN), KC = sum(KC)) %>% 
  group_by(geneId, taxid, pathId) %>% 
  filter(fullPath == max(fullPath)) %>%
  filter(accession %in% accession[rank == max(rank)][1]) %>%
  # filter(rank == max(rank)) %>% dplyr::slice(1) %>% 
  ungroup() 


allBact = allBact %>% 
  group_by(geneId, accession, taxid, genus, species, plasmid, LN, KC) %>% 
  summarise(extension = sum(pathScore[!start]) + bitScore[1],
            across(c(fullPath:maxOrder1), max), .groups = "drop") %>% 
  left_join(
    genesDetected %>%
      select(geneId, gene, subtype, cover, type), 
    by = "geneId"
  ) %>% distinct() %>% 
  mutate(pathScore = fullPath, depth = KC / LN) %>% 
  group_by(geneId) %>% 
  filter(extension >  0 | all(extension == 0)) %>% #Remove only start hits if any have extension
  ungroup()


#Adjust the scored for genes by presences of other genes in the same bact
# myData = genomeARG
# myData$taxid = myData$accession

adjustBact = function(myData, bactGroupStart = 0){
  
  if(nrow(myData) == 0) {
    return(list(myClusters = data.frame(), bactList = data.frame()))
  }
  
  if(n_distinct(myData$geneId) == 1 | n_distinct(myData$taxid) == 1){
    if(n_distinct(myData$taxid) == 1){
      
      #There is only one taxid
      bactList = myData %>% 
        select(taxid, genus, species, depth, pathScore) %>% 
        mutate(bactGroup = 1 + bactGroupStart, prob = 1, val = NA,
               runId, pipelineId = myPipelineId) %>% 
        filter(depth == max(depth)) %>% dplyr::slice(1) %>% 
        left_join(bactGenomeSize %>% select(genus, species, size) %>% distinct(),
                  by = c("genus", 'species')) 
      
      return(list(
        myClusters = data.frame(
          bactGroup = 1 + bactGroupStart,
          geneId = myData$geneId[1]),
        bactList = bactList))
    } else {
      
      #There is only one gene but multiple taxid
      geneMatrix = matrix(myData$pathScore / max(myData$pathScore), ncol = 1)
      colnames(geneMatrix) = unique(myData$geneId)
      rownames(geneMatrix) = myData$taxid
      myClusters = data.frame(
        taxid =  myData$taxid,
        bactGroup = 1 + bactGroupStart,
        primary = geneMatrix[,1] == 1,
        geneId = as.integer(unique(myData$geneId)),
        val = geneMatrix[,1])
    }
    
    
  } else {
    
    #There are multiple genes
    geneMatrix = myData %>%
      mutate(geneId = as.integer(geneId)) %>% 
      group_by(geneId, bact = taxid) %>%
      summarise(pathScore = round(max(pathScore),2), 
                extension = round(max(extension),2), 
                .groups = "drop") %>%
      left_join(genesDetected %>%
                  select(geneId, ARGgroup), by = "geneId")
    
    geneMatrix = geneMatrix %>%
      select(bact, geneId, pathScore) %>%
      pivot_wider(bact,
                  names_from = "geneId", values_from = "pathScore",
                  values_fill = 0) %>%
      column_to_rownames("bact") %>% as.matrix()
    
    #Normalise the matrix per column
    geneMatrix = apply(geneMatrix, 2, function(x) x / max(x))
    
    # Adjust by cover
    myAdjustment = data.frame(
      geneId = as.integer(colnames(geneMatrix))
    ) %>% left_join(
      genesDetected %>% select(geneId, cover), by = "geneId")
    
    #-----------
    # During adjustment make sure to weigh the influence of other ARG based
    # on the similarity in coverage. If the coverage is very different, the
    # weight of adjustment in lower (prevents FP association of low cover to high)
    ncols = ncol(geneMatrix)
    myAdjustment = (1 - abs(matrix(myAdjustment$cover, nrow = ncols, ncol = ncols) -
                              matrix(myAdjustment$cover, nrow = ncols, ncol = ncols, byrow = T)))
    
    
    myAdjustment = sapply(1:ncols, function(col){
      t(t(geneMatrix[,-col]) * myAdjustment[col,-col]) %>% rowSums() / 
        (ncols - 1)
    }) + 1
    
    #---------------
    ##Old adjustment techinique
    # geneMatrix = t(t(geneMatrix) * myAdjustment$cover)
    # 
    # #Calulcate the adjustment for each cell based on other scores in the row
    # # the higher the scores (e.g. other gened detected), the more scaled up
    # myAdjustment = (matrix(rowSums(geneMatrix), nrow = nrow(geneMatrix), 
    #                        ncol = ncol(geneMatrix)) - geneMatrix) /
    #   (ncol(geneMatrix) - 1) + 1
    
    #----------------
    
    geneMatrix = geneMatrix * myAdjustment
    geneMatrix = apply(geneMatrix, 2, function(x) x / max(x))
    
    
    myClusters = as_tibble(geneMatrix > 0, rownames = "taxid") %>%
      group_by(across(c(-taxid))) %>% 
      mutate(bactGroup = cur_group_id()) %>% ungroup()
    
    myClusters = as_tibble(geneMatrix, rownames = "taxid") %>% 
      left_join(myClusters %>% select(taxid, bactGroup), by = "taxid") 
    myClusters$primary = apply(myClusters[,!colnames(myClusters) %in% c("taxid", "bactGroup")], 1,
                               function(x) any(x == 1))
    
    #new ---
    x = myClusters %>% group_by(bactGroup) %>% filter(!all(!primary)) %>% 
      ungroup() 
    
    
    y = myClusters %>% group_by(bactGroup) %>% filter(all(!primary)) %>% 
      ungroup()
    # mutate(across(c(-taxid, -bactGroup, -primary), function(x) x > 0))
    
    
    a = apply(x %>% select(-taxid, -bactGroup, -primary), 1, function(z){
      
      apply(y %>% select(-taxid, -bactGroup, -primary), 1, function(a) cor(z,a))
      
    })
    
    colnames(a) = x$bactGroup
    
    y = cbind(a,y)
    y = y[,1:which(colnames(y) == "taxid")]
    
    y = y %>% pivot_longer(-taxid, names_to = "cluster") %>% 
      group_by(taxid) %>% filter(value == max(value)) %>% ungroup() %>% 
      mutate(cluster = str_remove(cluster, "\\.\\d+$"))
    
    myClusters = myClusters %>% left_join(y, by = "taxid") %>% 
      mutate(bactGroup = ifelse(is.na(cluster), bactGroup, cluster)) %>% 
      select(-cluster)
    #-----
    
    myClusters = myClusters %>% group_by(bactGroup) %>% filter(!all(!primary)) %>% 
      ungroup() %>% mutate(bactGroup = bactGroupStart + as.factor(bactGroup) %>% as.integer()) %>% 
      pivot_longer(c(-taxid, -bactGroup, -primary, -value), names_to = "geneId", values_to = "val") %>% 
      mutate(geneId = as.integer(geneId), taxid = as.integer(taxid)) %>% 
      filter(val > 0 | primary)
    
  }
  
  bactList = myClusters %>% group_by(bactGroup, taxid) %>% 
    summarise(prob = sum(val), .groups = "drop") %>% 
    group_by(bactGroup) %>% 
    mutate(prob = softmax(prob)) %>% ungroup() %>%
    mutate(taxid = as.integer(taxid), bactGroup = as.integer(bactGroup)) %>% 
    left_join(allBact %>% select(taxid, genus, species) %>% 
                distinct(), by = "taxid") %>% 
    group_by(bactGroup) %>% 
    filter(prob >= min(
      sort(unique(prob), decreasing = T)[1:(min(11, n_distinct(prob)))])
    ) %>% 
    mutate(val = prob / max(c(-Inf, prob), na.rm = T)) %>% 
    ungroup() %>% 
    arrange(bactGroup, desc(prob)) %>% ungroup() %>% 
    mutate(species = str_replace(species, "sp\\d+", "sp.")) %>% 
    left_join(bactGenomeSize %>% select(genus, species, size) %>% distinct(),
              by = c("genus", 'species')) %>% 
    mutate(runId = {{runId}},
           pipelineId = toProcess$pipelineId[sampleIndex]) %>% 
    left_join(
      myData %>% group_by(taxid, geneId) %>% 
        filter(fullPath == max(fullPath)) %>% dplyr::slice(1) %>% 
        group_by(taxid) %>% 
        summarise(
          depth = sum(KC)/sum(LN),
          pathScore = sum(pathScore), 
          .groups = "drop"
        ),
      by = "taxid"
    ) 
  
  return(list(myClusters = myClusters, bactList = bactList))
}


#Get all the ARG that are most likely found in genomes
# allBact$taxid = as.character(allBact$taxid)

genomeARG = allBact %>% 
  group_by(geneId, plasmid) %>% 
  summarise(path = max(fullPath), .groups = "drop") %>% 
  group_by(geneId) %>% 
  mutate(perc = path / max(path)) %>% 
  filter((plasmid & perc < 0.5) | !any(plasmid))

genomeARG = allBact %>% filter(geneId %in% genomeARG$geneId) #


#Check if there are any, proceed accordingly
if(nrow(genomeARG) > 0){
  
  #Adjust the bact presence across ARG
  genomeARG = adjustBact(genomeARG)
  
  #Get the clusters and bact list
  AMRclusters = genomeARG$myClusters %>%
    filter(val > 0) %>% 
    left_join(allBact %>% 
                select(taxid, geneId, pathScore), by = c("taxid", "geneId")) %>% 
    mutate(origin = "genome")
  
  AMRclusters = AMRclusters %>% 
    left_join(
      AMRclusters %>% 
        group_by(geneId, pathScore) %>% summarise(contenders = n(), .groups = "drop") %>% 
        group_by(geneId) %>% 
        arrange(desc(pathScore)) %>% 
        mutate(contenders = cumsum(contenders)) %>% ungroup(),
      by = c("geneId", "pathScore")
    )
  
  bactList = genomeARG$bactList %>% group_by(bactGroup) %>% 
    mutate(x = depth, depth = weighted.mean(depth, prob)) %>% ungroup()
}


#Check for likely plasmids
if(n_distinct(genomeARG$myClusters$geneId) < n_distinct(allBact$geneId)){
  
  #Check if any of the remaining ARG match a genome ARG species
  checkPlasmid = allBact %>% 
    filter(!geneId %in% genomeARG$myClusters$geneId) %>% 
    left_join(
      bactList %>% 
        # filter(val == 1) %>% #maybe remove?
        select(taxid, bactGroup, prob, genomeDepth = depth,
               genomePathscore = pathScore, val),
      by = "taxid"
    ) %>% left_join(
      AMRclusters %>% 
        # filter(val == 1) %>%
        select(taxid, geneId, primary) %>% 
        left_join(genesDetected %>% 
                    # mutate(cover = ifelse(type == "fragmentsOnly", cover2, cover)) %>% 
                    select(geneId, cover, startDepth, type),
                  by = "geneId") %>% 
        group_by(taxid) %>% 
        summarise(
          primary = all(primary),
          genomeCover = weighted.mean(cover, cover),
          genomeType = case_when(
            all(type == "noFragments") ~ "noFragments",
            all(type == "fragmentsOnly") ~ "fragmentsOnly",
            TRUE ~ "mixed"
          ), .groups = "drop"),
      by = "taxid"
    ) 
  
  checkPlasmid =  checkPlasmid %>% left_join(
    checkPlasmid%>% 
      group_by(geneId, pathScore) %>% summarise(contenders = n(), .groups = "drop") %>% 
      group_by(geneId) %>% 
      arrange(desc(pathScore)) %>% 
      mutate(contenders = cumsum(contenders)) %>% ungroup(),
    by = c("geneId", "pathScore")
  ) 
  
  #The closer in depth, the better
  # checkPlasmid = checkPlasmid %>% 
  #   mutate(bactGroup = ifelse(geneId == 143, NA, bactGroup)) %>% 
  #   group_by(geneId) %>% 
  #   mutate(x = max(0, fullPath[!is.na(bactGroup)]) >= 0.9*max(fullPath)) %>% 
  #   ungroup() %>% 
  #   filter((is.na(bactGroup) & !x) | (!is.na(bactGroup) & x)) %>% 
  #   select(-x) %>% rowwise() %>% 
  #   mutate(fullPath = fullPath / max(1, abs(depth - genomeDepth), na.rm = T)) 
  
  #Pathscore should be at in the top range
  checkPlasmid = checkPlasmid %>% 
    group_by(geneId) %>% 
    #Primary matches win if large enough
    mutate(x=any(primary & fullPath > 5000)) %>% 
    filter(fullPath <= ifelse(any(x), max(fullPath[primary & fullPath > 5000], na.rm = T), Inf)) %>% 
    mutate(
      #Give a slight edge to genome matches over plasmid
      fullPath = ifelse(!plasmid & !is.na(bactGroup), 
                        fullPath * 1.00, fullPath),
      #Pathscore should be at in the top range
      top = fullPath == max(fullPath)) %>% 
    ungroup() %>% select(-x)
  
  myPlasmids = list()
  temp = checkPlasmid
  #FInd the best matching plasmids
  while(nrow(temp > 0)){
    x = temp %>%
      filter(top) %>% 
      group_by(accession, genus, species) %>% 
      mutate(x = sum(fullPath)) %>% ungroup() %>% 
      filter(x == max(x))
    
    
    myPlasmids[[length(myPlasmids)+1]] = x
    
    
    temp = temp %>% filter(!geneId %in% x$geneId) 
  }
  
  #Link genome taxid to the plasmids and asssign if match criteria
  # sapply(myPlasmids, function(x) "2518" %in%  x$geneId)
  # x = myPlasmids[[2]]
  temp = map_df(myPlasmids, function(x){
    
    test = checkPlasmid %>% filter(geneId %in% x$geneId) %>% 
      group_by(taxid, geneId) %>% filter(fullPath == max(fullPath)) %>% 
      dplyr:: slice(1) %>% 
      group_by(taxid) %>% mutate(x = sum(fullPath)) %>% 
      ungroup() %>% mutate(x = x / max(x)) %>% 
      group_by(taxid, geneId) %>% 
      summarise(fullPath = sum(fullPath), depth = mean(depth),
                across(c(bactGroup:genomeType), max), .groups = "drop") %>%
      mutate(fullPath = ifelse(is.na(bactGroup), fullPath, fullPath*1.05)) %>% 
      group_by(geneId) %>% 
      filter(depth >= 0.75 * genomeDepth | is.na(genomeDepth) |
               abs(depth - genomeDepth) < 10) %>%
      #In case top is not in bactgroup, but first one in one has total score > 5000,
      #Remove the top matches until the first one with bactgroup
      filter(fullPath <= ifelse(max(0, fullPath[!is.na(bactGroup)]) > 5000,
                                max(fullPath[!is.na(bactGroup)]), max(fullPath))) %>% 
      filter(fullPath == max(fullPath))

    
  }, .id = "plasmidGroup") %>% 
    group_by(plasmidGroup) %>% 
    filter((bactGroup %in% bactGroup[!is.na(bactGroup)])| 
             all(is.na(bactGroup)))
  
  #Build data frame with the annotated plasmids
  myPlasmids = map_df(myPlasmids, function(x) x, .id = "plasmidGroup") %>% 
    select(plasmidGroup, geneId, bestAccession = accession) %>% 
    left_join(temp %>% select(plasmidGroup, taxid), by = "plasmidGroup") %>% 
    left_join(checkPlasmid, by = c("geneId", "taxid")) 
  
} else {
  
  myPlasmids = data.frame(
    accession = character(), geneId = integer(), bactGroup = integer(), 
    plasmidGroup = integer(), pathScore = numeric())
}

#Certainty can also be derived from complexity of gfa, the more segments and links
#the more complex and thus the more possible confusion
# score 1 = two paths, best scenaria
# score >1 = one path, actually fragment ...
# score < 1, the smaller the more paths, the more complex
pathComplexity = pathData %>% group_by(geneId) %>% 
  summarise(pathCompl = round(2/n_distinct(pathId), 3)) %>% 
  mutate(geneId = as.integer(geneId))

#Link plasmids to genomes if possible
test = bind_rows(
  allBact %>% 
    left_join(
      AMRclusters %>% filter(val == 1) %>% 
        select(taxid, geneId, bactGroup),
      by = c("taxid", "geneId")
    ) %>% filter(!is.na(bactGroup)) %>% 
    group_by(geneId, taxid) %>% filter(fullPath == max(fullPath)) %>% 
    dplyr::slice(1) %>% ungroup() %>% distinct(),
  
  #Plasmid
  allBact %>% 
    left_join(
      myPlasmids %>% distinct() %>% 
        # group_by(plasmidGroup, geneId) %>% 
        # filter(pathScore >= 0.75*max(pathScore, na.rm = T)) %>% 
        # filter(all(!bactGroup %in% genomeARG$myClusters$bactGroup) | 
        #          bactGroup %in% genomeARG$myClusters$bactGroup) %>% 
        group_by(bactGroup, geneId) %>% 
        filter(pathScore == max(0, pathScore)) %>% dplyr::slice(1) %>% 
        select(accession, geneId, bactGroup, plasmidGroup) %>% distinct(),
      by = c("accession", "geneId")
    ) %>% filter(!is.na(plasmidGroup)) 
  
) %>% left_join(pathComplexity, by = "geneId") %>% 
  select(bactGroup, genus, species, geneId, gene, subtype, type, pathCompl, 
         fullPath, everything()) %>% 
  arrange(bactGroup, gene, subtype) %>% mutate(pipelineId = myId)

# return(test)

# }, error = function(x){
#   return(data.frame(pipelineId = myId))
# })
# }

sample

