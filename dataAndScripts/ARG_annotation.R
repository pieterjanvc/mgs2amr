#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse, warn.conflicts = FALSE)
library(RSQLite)

maxPathDist = 1500

#Get the arguments
tempFolder = formatPath("temp/", endWithSlash = F)
pipelineIds = NULL

#Load the ARG and the sample list to process
myConn = dbConnect(SQLite(), "dataAndScripts/meta2amr.db")

ARG = dbReadTable(myConn, "ARG") %>% 
  mutate(geneId = as.character(geneId))

toProcess = dbReadTable(myConn, "pipeline") %>% 
  filter(statusCode == 4)
if(length(pipelineIds) > 0){
  toProcess = toProcess %>% filter(pipelineId %in% pipelineIds)
}

dbDisconnect(myConn)


# ---- FUNCTIONS ----
#********************
cutOff = function(numbers, percent = 0.95){
  if(length(numbers) < 2){
    return(numbers)
  } else if(length(numbers) == 2){
    return(numbers[numbers >= max(numbers) * percent])
  }
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique() %>% min()
}

sample = "temp/E12Heidi011_SRR4017917_0.1_1615998709"
sample = "temp/test1"
sample = toProcess$tempFolder[1]

#Get the isolate info
# myConn = dbConnect(SQLite(), "sideStuff/isolateBrowserData.db")
# 
# sampleInfo = str_extract(sample, "SRR[^\\_]+")
# 
# sampleInfo = dbReadTable(myConn, "sampleInfo") %>% filter(Run == sampleInfo) %>% 
#   select(biosample_acc, Run, bact) %>% 
#   left_join(dbReadTable(myConn, "genoTypes"), by = "biosample_acc")
# 
# dbDisconnect(myConn)
options(readr.num_columns = 0)
i = 6
results = list()
results = map_df(1:nrow(toProcess), function(i){
# for(i in 1:2){
  
  sample = toProcess$tempFolder[i]
  sampleName = str_extract(sample, "[^\\/]+(?=_\\d+$)")
  genesDetected = read_csv(paste0(sample, "/genesDetected/genesDetected.csv")) %>% 
    mutate(subtype = as.character(subtype))
  
  print(paste("Annotating the genes for", sampleName))
  
  # ---- FILTERING DATA ----
  #*************************
  
  #Read all blast output
  blastOut = lapply(
    list.files(sample, full.names = T,
               pattern = "blastSegmentsClustered\\d+.json.gz"),
    blast_readResults, outFormat = "dataFrame1") %>% bind_rows()
  
  #Extract data we need + transform
  blastOut = blastOut %>% 
    select(query_title, hitId, taxid, accession , bact = title, bit_score, 
           score, evalue, identity, query_len, query_from, 
           query_to, hit_from, hit_to, align_len) %>% 
    mutate(bact = str_remove_all(bact, "[^\\w\\s]")) %>% 
    extract(bact, c("genus", "species", "extra"), 
            regex = "(\\w+)\\s+(\\w+)($|\\s+.*)") %>% 
    #Sp. will be pasted with taxid to make it unique
    mutate(species = ifelse(species == "sp", paste0(species, taxid), species)) %>% 
    filter(!species %in% c("bacterium", "xxx", "sp") & 
             !is.na(species)) %>% 
    mutate(
      subspecies = str_extract(extra, "(?<=subsp\\s)[^\\s]+"),
      strain = str_extract(extra, "(?<=strain\\s)[^\\s]+") %>% 
        str_remove(" chromosome.*| plasmid.*| complete.*"),
      plasmid = str_detect(extra, "plasmid|Plasmid"),
      plasmidName = str_extract(extra, "(?<=plasmid\\s)[^\\s]+")
    )

  #Expand results from clustering segments before blast
  clusterOut = read.table(paste0(sample, "/blastSegments.out")) %>% 
    select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
    mutate(
      clusterId = ifelse(clusterId == "*", segmentId, clusterId),
      geneId = str_extract(segmentId, "^([^_]+)")
      ) %>% 
    filter(clusterId %in% blastOut$query_title) 
  
  blastOut = clusterOut %>% 
    left_join(blastOut, by = c("clusterId" = "query_title")) %>% 
    mutate(start = str_detect(segmentId, "_start$"))
  
  
  #Get the KC and LN for the segments
  segmentInfo =  read_csv(paste0(sample, "/blastSegments.csv"), col_types = cols()) %>% 
    select(-sequence, -name, -geneId)
  
  #Check the bacteria
  allBact = blastOut %>% 
    select(-c(score:hit_to)) %>% 
    # filter(!plasmid, species != "sp") %>% 
    left_join(segmentInfo, by = c("segmentId" = "blastId")) %>% #Add segment info
    group_by(start, segmentId) %>% 
    filter(bit_score == max(bit_score)) %>% 
    group_by(start, segmentId, plasmid) %>%
    mutate(maxBit = max(bit_score)) %>% 
    group_by(start, segmentId) %>% 
    filter(maxBit == max(maxBit), !plasmid) %>% 
    mutate(onlyOne = length(unique(genus)) == 1 & 
             length(unique(species))) %>% 
    filter(onlyOne) %>% 
    group_by(start, genus, species) %>% 
    summarise(
      bit_score = sum(bit_score),
      depth = sum(KC / sum(LN)),
      .groups = "drop"
    ) %>% 
    group_by(start) %>% 
    filter(depth <= depth[bit_score == max(bit_score)])
  
  
  pathData = list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T)
  
  pathData = map_df(pathData, function(myGFA){
    gfa = gfa_read(myGFA)
    if(nrow(gfa$links) > 0){
      segmentOfInterest = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
        filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% slice(1) %>% pull(name)
      
      x = gfa_pathsToSegment(gfa, segmentOfInterest, returnList = T, pathSegmentsOnly = T) %>% 
        map_df(function(path){
          data.frame(
            pathId = path$id,
            startOrientation = path$startOrientation,
            segmentId = path$segmentOrder)
        })
      
      if(nrow(x) == 0 ){
        data.frame()
      }
      else{
        x %>% 
          left_join(
            pathsToSegmentTable(gfa, segmentOfInterest) %>% 
              filter(dist < Inf) %>% group_by(segment) %>% 
              filter(dist == min(dist)) %>% 
              select(segment, KC, LN, dist) %>% 
              mutate(geneId = str_extract(segment, "^\\d+")) %>% 
              distinct(),
            by = c("segmentId" = "segment")
          ) %>% 
          filter(LN >=250) %>% 
          group_by(geneId, pathId) %>% mutate(
            pathId = as.integer(pathId),
            order = n():1,
            type = "full")
      }
      
      
      
    } else {
      data.frame()
    }
    
  })
  
  
  fragments = gfa_read(paste0(sample, "/fragmentsOnly.gfa"))
  # fragments = gfa_read(paste0(sample, "/fragmentGFA.gfa"))
  if(nrow(fragments$segments) > 0){

    fragments = fragments$segments %>% 
      select(segmentId = name, LN, KC, geneId) %>% 
      group_by(geneId) %>% 
      mutate(pathId = 1, startOrientation = 0:(n()-1),
             order = 1, type = "fragment", dist = 0)
  } else {
    fragments = data.frame()
  }
  
  pathData = bind_rows(pathData, fragments)
  
  #Get the result
  result = 
    blastOut %>% 
    # mutate(species = ifelse(species == "cloacae", "hormaechei", species)) %>% 
    filter(
      extra != "",
      genus %in% allBact$genus, #Only keep the bact of interest
      species %in% allBact$species
    ) %>% #fix weird extra
    select(-hit_from, -hit_to) %>% #somtimes more than one match like repeat? (ignore)
    mutate(
      coverage = align_len / query_len, #coverage of segment
      extra = ifelse(species == "sp", "xxx", extra) #fix species things to be more general
    ) %>%
    filter(coverage >= 0.95 | (align_len > 500 & coverage >= 0.90)) %>% #filter by segment coverage
    group_by(start, segmentId, geneId, genus, species, extra, plasmid) %>% 
    filter(bit_score == max(bit_score)) %>% distinct() %>% #only one subspecies per segment
    ungroup() %>% left_join(pathData, by = c("geneId", "segmentId")) %>% 
    filter(!is.na(dist)) %>% rowwise() %>% #Only keep segments that are in a path to start
    mutate(
      dist = ifelse(dist < 1, 1, dist) #Avoid division by 0
      ) %>% group_by(geneId, accession, genus, species, extra, pathId, 
                     start, startOrientation) %>% 
    mutate(
      completePath = length(unique(order[order != 1])) == (max(order) - 1),
      score = sum(score),
      bit_score = sum(bit_score)
      ) %>%
    ungroup() %>% 
    filter(completePath) %>% select(-dist, -order, -hitId) %>% 
    distinct()
  
  result$startOrientation[result$start] = 0
  result = result %>% 
    select(start,geneId, accession, genus, species, 
           extra, plasmid, pathId, bit_score, score) %>% 
    group_by(geneId, genus, species, extra, start, pathId) %>% 
    mutate(
      score = max(score),
      bit_score = max(bit_score)
    ) %>% 
    arrange(desc(bit_score)) %>% slice(1) %>% 
    group_by(geneId, accession, genus, species, extra, pathId) %>% 
    mutate(
      total_score = sum(score),
      total_bit_score = sum(bit_score)
    ) %>% ungroup() %>% 
     left_join(ARG %>% select(geneId, gene, subtype, clusterNr), by = "geneId") %>% #Add ARG info
    mutate(val = total_bit_score)
  
    
  #Aggregate on species level
  test = result %>% 
    group_by(geneId, gene, subtype, accession, genus, species, extra, plasmid, pathId) %>% 
    summarise(startVal = max(0,bit_score[start]), 
              notStartVal = max(0,bit_score[!start]),
              val = total_bit_score[1]) %>% 
    group_by(geneId, gene, subtype, genus, species) %>% 
    mutate(
      startVal = ifelse(startVal == 0, max(startVal) / 2, startVal),
      notStartVal = ifelse(notStartVal == 0, max(notStartVal) / 2, notStartVal),
      val = startVal + notStartVal
    ) %>% 
    group_by(geneId, pathId) %>%
    filter(val == max(val)) %>% 
    group_by(geneId) %>%
    filter(val == max(val)) %>% ungroup() %>% 
    select(-pathId) %>% distinct()
  
  details = test
  test = test %>% 
    group_by(geneId, gene, subtype, genus, species, plasmid) %>% 
    arrange(desc(val)) %>% slice(1) %>% 
    group_by(geneId, gene, subtype, genus, species) %>%
    mutate(plasmid = case_when(
      sum(any(plasmid) + any(!plasmid)) == 2 ~ max(val[plasmid]) / 
        (max(val[plasmid]) + max(val[!plasmid])),
      any(plasmid) ~ 1,
      TRUE ~ 0
    ))
  
  details = result
  
  test = test %>% group_by(geneId, gene, subtype, genus, species) %>% 
    arrange(desc(val)) %>% slice(1) %>% 
    ungroup() %>% select(-extra, -accession) %>% 
    left_join(genesDetected %>% 
                mutate(geneId = as.character(geneId)) %>% 
                select(geneId, nBases, startPerc, startDepth, cover1, 
                       cover2, KCsum, type),
              by = "geneId") %>% 
    rowwise() %>% 
    mutate(
      pipelineId = toProcess$pipelineId[i],
      startVal = min(startVal / (nBases*1.81), 1),
      notStartVal = min(notStartVal / (maxPathDist*2*1.81), 1)
    ) %>% 
    arrange(genus, species, desc(val))
  
  test
  
  # write_csv(result, paste0(sample, "/annotation.csv"))
})

results = results %>% select(pipelineId, everything())
write_csv(results, 'sideStuff/annotationResults_2.csv')

