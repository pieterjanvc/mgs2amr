#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse)
library(RSQLite)

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
cutOff = function(numbers){
  if(length(numbers) < 2){
    return(numbers)
  }
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique() %>% min()
}

sample = "temp/G2Heidi019_1612974558"
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
  
  #Extract data we need 
  blastOut = blastOut %>% 
    select(query_title, hitId, taxid, accession , bact = title, bit_score, 
           score, evalue, identity, query_len, query_from, 
           query_to, hit_from, hit_to, align_len) %>% 
    mutate(bact = str_remove_all(bact, "[^\\w\\s]")) %>% 
    extract(bact, c("genus", "species", "extra"), 
            regex = "(\\w+)\\s+(\\w+)($|\\s+.*)") %>% 
    filter(!species %in% c("bacterium", "xxx") & !is.na(species)) %>% 
    mutate(plasmid = str_detect(extra, "plasmid|Plasmid"))
  
  #Expand results from clustering segments before blast
  clusterOut = read.table(paste0(sample, "/blastSegments.out")) %>% 
    select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
    mutate(clusterId = ifelse(clusterId == "*", segmentId, clusterId)) %>% 
    extract(segmentId, into = c("geneId", "segment"), regex = "^([^_]+)_(.*)", remove = F) %>% 
    filter(clusterId %in% blastOut$query_title) %>% arrange(clusterId)
  
  blastOut = clusterOut %>% 
    left_join(blastOut, by = c("clusterId" = "query_title"))
  
  
  #Get the KC and LN for the segments
  segmentInfo =  read_csv(paste0(sample, "/blastSegments.csv")) %>% 
    select(-sequence, -name, -geneId, -CL)
  
  #Check the bacteria
  allBact = blastOut %>% select(-c(score:hit_to)) %>% 
    # filter(!plasmid, species != "sp") %>% 
    left_join(segmentInfo, by = c("segmentId" = "blastId")) %>% #Add segment info
    group_by(segmentId) %>% 
    filter(bit_score == max(bit_score)) %>% 
    group_by(segmentId, plasmid) %>%
    mutate(maxBit = max(bit_score)) %>% 
    group_by(segmentId) %>% 
    filter(maxBit == max(maxBit), !plasmid) %>% 
    mutate(onlyOne = length(unique(genus)) == 1 & 
             length(unique(species))) %>% 
    filter(onlyOne) %>% 
    group_by(genus, species) %>% 
    summarise(
      bit_score = sum(bit_score),
      depth = sum(KC / sum(LN)),
      .groups = "drop"
    ) %>% 
    filter(depth <= depth[bit_score == max(bit_score)])
  
  #Get all data from the paths
  test = list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T)
  
  pathData = map_df(test, function(myGFA){
    gfa = gfa_read(myGFA)
    segmentOfInterest = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
      filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% slice(1) %>% pull(name)
    
    if(nrow(gfa$links) > 0){
      pathsToSegmentTable(gfa, segmentOfInterest) %>% 
        filter(dist < Inf) %>% group_by(segment) %>% 
        filter(dist == min(dist)) %>% 
        select(segment, KC, dist) %>% 
        mutate(geneId = str_extract(myGFA, "\\d+(?=_simplified.gfa)"))
    } else {
      data.frame()
    }
    
  })
  
  #Get the result
  blastOut = 
    blastOut %>% 
    filter(
      extra != "",
      genus %in% allBact$genus,
      species %in% allBact$species
    ) %>% #fix weird extra
    select(-hit_from, -hit_to) %>% #somtimes more than one match like repeat? (ignore)
    mutate(
      coverage = align_len / query_len, #coverage of segment
      extra = ifelse(species == "sp", "xxx", extra) #fix species things to be more general
    ) %>%
    filter(coverage >= 0.9 | align_len > 250) %>%
    group_by(segmentId, geneId, genus, species, extra, plasmid) %>% 
    filter(bit_score == max(bit_score)) %>% distinct() %>% 
    ungroup() %>% left_join(pathData, by = c("geneId", "segment")) %>% 
    filter(!is.na(dist)) %>% rowwise() %>% 
    mutate(
      dist = dist + 1,
      val = sum(1 / dist:(dist + align_len)) * bit_score / align_len) %>% 
    group_by(geneId, accession, genus, species, extra, plasmid) %>% 
    summarise(val = sum(val)) %>% 
    group_by(geneId, genus, species, plasmid) %>% 
    filter(val == max(val)) %>% ungroup() %>% 
    left_join(ARG %>% select(geneId, gene, subtype, clusterNr), by = "geneId") %>% 
    select(gene, subtype, clusterNr, accession, genus, species, extra, plasmid, val) %>% distinct() %>% 
    group_by(gene, subtype) %>%  filter(val >= cutOff(val))

    write_csv(blastOut, paste0(sample, "/annotation.csv"))
}
