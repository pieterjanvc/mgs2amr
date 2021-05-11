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
sample = "temp/E12Heidi011_SRR4017912_0.05_1615998142"
sample = toProcess$tempFolder[12]

#Get the isolate info
myConn = dbConnect(SQLite(), "sideStuff/isolateBrowserData.db")

sampleInfo = str_extract(sample, "SRR[^\\_]+")

sampleInfo = dbReadTable(myConn, "sampleInfo") %>% filter(Run == sampleInfo) %>% 
  select(biosample_acc, Run, bact) %>% 
  left_join(dbReadTable(myConn, "genoTypes"), by = "biosample_acc")

dbDisconnect(myConn)

genesDetected = read_csv(paste0(sample, "/genesDetected/genesDetected.csv")) %>% 
  mutate(subtype = as.character(subtype))

for(sample in toProcess$tempFolder){
  
  sampleName = str_extract(sample, "[^\\/]+(?=_\\d+$)")
  print(paste("Annotating the genes for", sampleName))
  
  # ---- FILTERING DATA ----
  #*************************
  
  #Read all blast output
  blastOut = lapply(
    list.files(sample, full.names = T,
               pattern = "blastSegmentsClustered\\d+.json.gz"),
    blast_readResults, outFormat = "dataFrame1") %>% bind_rows()
  
  # test = blastOut #REMOVE !!!
  # blastOut = test
  
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
    mutate(clusterId = ifelse(clusterId == "*", segmentId, clusterId)) %>% 
    extract(segmentId, into = c("geneId", "segment"), regex = "^([^_]+)_(.*)", remove = F) %>% 
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
    filter(depth <= depth[bit_score == max(bit_score)])
  
  
  pathData = list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T)
  
  pathData = map_df(pathData, function(myGFA){
    gfa = gfa_read(myGFA)
    if(nrow(gfa$links) > 0){
      segmentOfInterest = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
        filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% slice(1) %>% pull(name)
      
      gfa_pathsToSegment(gfa, segmentOfInterest, returnList = T, pathSegmentsOnly = T) %>% 
        map_df(function(path){
          data.frame(
            pathId = path$id,
            startOrientation = path$startOrientation,
            segment = path$segmentOrder)
        }) %>% 
        left_join(
          pathsToSegmentTable(gfa, segmentOfInterest) %>% 
            filter(dist < Inf) %>% group_by(segment) %>% 
            filter(dist == min(dist)) %>% 
            select(segment, KC, LN, dist) %>% 
            mutate(geneId = str_extract(segment, "^\\d+")) %>% 
            distinct(),
          by = "segment"
        ) %>% 
        filter(LN >=250) %>% 
        group_by(geneId, pathId) %>% mutate(
          pathId = as.integer(pathId),
          order = n():1,
          type = "full")
      
      
    } else {
      data.frame()
    }
    
  })
  
  fragments = gfa_read(paste0(sample, "/fragmentGFA.gfa"))
  fragments = fragments$links %>% select(from, to, geneId) %>% 
    mutate(pathId = 1:n()) %>% 
    pivot_longer(c(from, to), names_to = NULL, values_to = "segment") %>% 
    left_join(fragments$segments %>% select(name, LN, KC), by = c("segment" = "name")) %>% 
    group_by(pathId) %>% 
    filter(any(LN >= 250)) %>% ungroup() %>% 
    mutate(
      dist = ifelse(str_detect(segment, "_start$"), -1* LN, 0),
      order = ifelse(dist == 0, 1, 2),
      type = "fragment") 
  
  fragments = fragments %>% left_join(
    fragments %>% select(pathId, segment) %>% 
      filter(str_detect(segment, "_start$")) %>% 
      group_by(segment) %>% mutate(startOrientation = 1:n() -1),
    by = c("pathId", "segment")) %>% group_by(pathId) %>% 
    mutate(startOrientation = sum(startOrientation, na.rm = T))
  
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
    ungroup() %>% left_join(pathData, by = c("geneId", "segment")) %>% 
    filter(!is.na(dist)) %>% rowwise() %>% #Only keep segments that are in a path to start
    mutate(
      dist = ifelse(dist < 1, 1, dist) #Avoid division by 0
      # val = bit_score
        # sum(1 / dist:(dist + align_len)) * 
        # bit_score / align_len #Scale the bit_score based off distance to start
      # val2 = sum(1 / dist:(dist + min(align_len, 1500))) * 
      #   min(score, 1500) / min(align_len, 1500)
      ) %>% group_by(geneId, accession, genus, species, extra, pathId, 
                     start, startOrientation) %>% 
    # mutate(completePath = sum(unique(order)) / sum(min(2,max(order)):max(order))) %>% 
    # mutate(completePath = sum(unique(order)) / sum(1:max(order))) %>% 
    mutate(
      # completePath = length(unique(order[order != 1])) == (max(order) - 1),
      completePath = length(unique(order[order != 1])) == (max(order) - 1),
      score = sum(score),
      bit_score = sum(bit_score)
      ) %>%
    ungroup() %>% 
    filter(completePath) %>% select(-dist, -order, -hitId) %>% 
    distinct()
  
  # test = result %>% filter(!start) %>% 
  #   select(segmentId, genus, species,val) %>% 
  #   group_by(segmentId, genus, species) %>% 
  #   filter(val == max(val)) %>% 
  #   distinct() %>% group_by(segmentId) %>% 
  #   filter(val >= cutOff(val, 0.99)) %>% 
  #   group_by(segmentId) %>% 
  #   mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% 
  #   mutate(geneId = as.integer(str_extract(segmentId, "^\\d+"))) %>% 
  #   group_by(geneId) %>% 
  #   filter(length(unique(paste(genus, species))) == 1)
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
    # group_by(geneId, pathId) %>% 
    # filter(total_bit_score == max(total_bit_score))
    # summarise(val = sum(val), .groups = "drop") %>% #Summary per geneId
    # group_by(start,geneId, genus, species, plasmid) %>% 
    # filter(total_bit_score == max(total_bit_score)) %>% ungroup() %>% #only keep best subspecies (extra) per species
    left_join(ARG %>% select(geneId, gene, subtype, clusterNr), by = "geneId") %>% #Add ARG info
    # select(geneId, gene, subtype, clusterNr, accession, genus, species, extra, 
    #        start, plasmid, bit_score, score, total_bit_score) %>% 
    # distinct() %>% 
    mutate(val = total_bit_score)
  
    #Keep the details for subspecies separate
  # details = result %>% 
  #   select(gene, subtype, accession, genus, species, extra, start, plasmid, val)
    
    #Aggregate on species level
  test = result %>% 
    # group_by(gene, subtype, genus, species, extra, start, plasmid) %>% 
    # arrange(desc(val)) %>% slice(1) %>% select(-accession, -extra) %>% 
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
    # group_by(geneId, gene, subtype, accession, genus, species, extra, plasmid) %>% 
    # mutate(
    #   startVal = startVal[1],
    #   notStartVal = sum(notStartVal),
    #   ) %>% 
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
      startVal = min(startVal / (nBases*1.81), 1),
      notStartVal = min(notStartVal / (maxPathDist*2*1.81), 1)
    ) %>% 
    arrange(genus, species, desc(val))
  
  # test = result %>% select(-accession, -extra) %>% distinct() %>% 
  #     group_by(start,gene, subtype) %>%  filter(val >= cutOff(val)) %>%
  #     # group_by(start,gene, subtype, genus, species) %>%  filter(val >= cutOff(val)) %>% 
  #     left_join(genesDetected, by = c("gene", "subtype")) %>% 
  #     left_join(test %>% select(-segmentId, -val, filter = n) %>% distinct(), 
  #               by = c("geneId", "genus", "species")) %>% 
  #     group_by(geneId) %>% 
  #     mutate(
  #       filter = ifelse(is.na(filter) & any(filter == 1,  na.rm = T), 0, 1),
  #       cover = ifelse(type == "fragmentsOnly", cover2, cover1)) %>% 
  #     filter(filter == 1) %>% 
  #     group_by(geneId, genus, species) %>% filter(val == max(val)) %>% ungroup() %>% 
  #     group_by(gene) %>% filter(cover >= min(0.9, max(cover)) | type != "fragmentsOnly")
  #   
  # result = result %>% 
  #     select(geneId, gene, subtype, genus, species, plasmid, val, val2, cover, startPerc, type, everything()) %>%  
  #     arrange(gene, desc(cover), desc(val), desc(startPerc),genus, species)
  #   
  # result = result %>% 
  #     select(geneId, gene, subtype, genus, species, plasmid, val, val2, cover, startPerc, type, everything()) %>%  
  #     arrange(genus, species, desc(type), gene, desc(cover), desc(val), desc(startPerc))
    
  
  write_csv(result, paste0(sample, "/annotation.csv"))
}



