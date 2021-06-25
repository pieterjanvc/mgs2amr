#**************************
# ---- ARG ANNOTATION ----
#*************************

library(gfaTools)
library(tidyverse, warn.conflicts = FALSE)
library(RSQLite)
library(igraph)
library(visNetwork)

#Get the arguments
baseFolder = "/mnt/meta2amrData/meta2amr/"
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

#Generate a list out of the settings file
settings = readLines(paste0(baseFolder, "settings.txt"))
settings = settings[str_detect(settings,"^\\s*[^=#]+=.*$")]
settings = str_match(settings, "\\s*([^=\\s]+)\\s*=\\s*(.*)")
settings = setNames(str_trim(settings[,3]), settings[,2])


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
i = 30
results = list()
# results = map_df(1:nrow(toProcess), function(i){
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
    filter(!genus %in% c("uncultured", "mixed"),
           !species %in% c("bacterium") & 
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
  allBact2 = blastOut %>% 
    rowwise() %>% 
    mutate(
      ident = identity / align_len,
      coverage = min(align_len / query_len, 1),
    ) %>%
    filter(coverage == 1 , ident >= 0.95) %>%
    select(segmentId, geneId, hitId, bit_score, 
           accession, taxid, genus, species, extra) %>% 
    group_by(segmentId, geneId, hitId, accession) %>% 
    filter(bit_score == max(bit_score)) %>% #Remove multiple hits (keep best)
    slice(1) %>% ungroup() %>% 
    
    left_join(segmentInfo, by = c("segmentId" = "blastId")) %>% #Add segment info
    
    group_by(segmentId) %>% 
    mutate(onlyOne = length(unique(genus)) < 2) %>% 
    filter(onlyOne) %>% 
    group_by(taxid, accession, genus, species, extra) %>% 
    summarise(
      bit_score = sum(bit_score),
      depth = sum(KC) / sum(LN),
      .groups = "drop"
    ) %>% 
    group_by(taxid, genus, species) %>% 
    filter(bit_score == max(bit_score)) %>% 
    slice(1) %>% ungroup()
  
  blastOut = blastOut %>% 
    filter(taxid %in% allBact2$taxid)
  
  
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

  if(nrow(fragments$segments) > 0){
    
    fragments = fragments$segments %>% 
      select(segmentId = name, LN, KC, geneId) %>% 
      group_by(geneId) %>% 
      mutate(pathId = 1, startOrientation = 0:(n()-1),
             order = 1, type = "fragment", dist = -1) #change back !!!?
  } else {
    fragments = data.frame()
  }
  
  pathData = bind_rows(pathData, fragments)  %>%  ungroup()
  

  pathData = pathData %>% group_by(geneId, pathId, startOrientation) %>% 
    mutate(
      depth = KC / LN,
      totalLN = sum(LN[dist >= 0])) %>% ungroup()
  
  
  result = blastOut %>% 
    rowwise() %>% mutate(
      ident = identity / align_len,
      coverage = max(align_len / query_len, 1),
    ) %>%
    filter(coverage >= 0.9, ident >= 0.95) %>% 
    group_by(segmentId, geneId, hitId, accession) %>% 
    filter(bit_score == max(bit_score)) %>% #Remove multiple hits (keep best)
    slice(1) %>% ungroup() %>% 
    select(segmentId, geneId, accession, taxid, plasmid, genus, species, strain, 
           extra, bit_score, coverage, ident) %>% 
    left_join(pathData %>% select(segmentId, pathId, order, startOrientation, 
                                  dist, depth, type) %>% 
                distinct(), by = "segmentId") %>% 
    filter(!is.na(dist), order > 0) %>% 
    left_join(genesDetected %>% 
                mutate(geneId = as.character(geneId)) %>% 
                select(geneId, gene, subtype), by = "geneId") %>% 
    rowwise() %>% 
    mutate(pathScore = bit_score * coverage / (max(dist,0) * 0.01 + 1))
  
  result = result %>% 
    group_by(accession, pathId) %>% 
    mutate(n = n()) %>% 
    group_by(geneId, accession, pathId) %>% 
    group_by(geneId, gene, subtype, accession, taxid, genus, species, 
             plasmid, pathId, type) %>% 
    summarise(
      pathScore = sum(pathScore),
      depth = mean(depth), .groups = "drop") %>% 
    group_by(geneId) %>% filter(pathScore == max(pathScore)) %>% 
    select(-pathId) %>% distinct()
  
  #---- Remove duplicate genes ---
  #-------------------------------
  getSegments = map_df(list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T), function(x) gfa_read(x)$segments)
  
  myFilter = pathData %>% filter(order < 3) %>% pull(segmentId) %>% unique()
  
  getSegments = getSegments %>% filter(name %in% myFilter)
  
  fasta_write(getSegments$sequence, 
              sprintf("%s/ARGsim.fasta", sample),
              getSegments$name, type = "n")
  
  #Add the reverse complement (usearch does not do that when calc_distmx)
  system(sprintf("%s -fastx_revcomp %s -label_suffix _RC -fastaout %s >/dev/null 2>&1",
                 settings["usearch"],
                 sprintf("%s/ARGsim.fasta", sample),
                 sprintf("%s/ARGsim_RC.fasta", sample)))
  system(sprintf("cat %1$s %2$s > %3$s; rm %1$s %2$s >/dev/null 2>&1",
                 sprintf("%s/ARGsim.fasta", sample),
                 sprintf("%s/ARGsim_RC.fasta", sample),
                 sprintf("%s/ARGsimilarities.fasta", sample)))
  
  #Use cluster_fast to reduce number of segments by grouping in identity clusters
  verbose = 2
  system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.1 %s",
                 settings["usearch"],
                 sprintf("%s/ARGsimilarities.fasta", sample),
                 sprintf("%s/ARGsimilarities.out", sample),
                 ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))
  
  
  #Group identical ARG and remove the duplicate (weak) ones
  identicalARG = read_delim(sprintf("%s/ARGsimilarities.out", sample), "\t", col_names = F) %>% 
    filter(X1 != X2) %>% 
    mutate(X1 = str_remove(X1, "_RC$"), X2 = str_remove(X2, "_RC$")) %>% 
    left_join(pathData %>% filter(order < 3) %>% 
                select(segmentId, d1 = depth, o1 = order, LN1 = LN, KC1 = KC) %>% 
                distinct(), by = c("X1" = "segmentId")) %>% 
    left_join(pathData %>% filter(order < 3) %>% 
                select(segmentId, d2 = depth, o2 = order, LN2 = LN, KC2 = KC) %>% 
                distinct(), by = c("X2" = "segmentId")) %>% 
    filter(o1 == o2) %>% 
    rowwise() %>% 
    mutate(
      gene1 = str_extract(X1, "^\\d+"),
      gene2 = str_extract(X2, "^\\d+"),
      LNdiff = min(LN1, LN2) / max(LN1, LN2)
    ) %>% 
    #Remove alignments that are too big in size difference, unless start
    filter(LNdiff > 0.9 | (X3 == 0 & o1 == 1)) %>% 
    mutate(gene = ifelse(d1 > d2, gene1, gene2),
           depth = ifelse(d1 > d2, d1, d2))
  
  #Build a graph to see how the genes are connected
  myFilter = graph_from_data_frame(data.frame(
    from = identicalARG$gene1, to = identicalARG$gene2
  ), directed = F)
  
  #Only consider similarities if at least two pieces align
  myFilter = as_adjacency_matrix(myFilter, sparse = F)
  myFilter = as.data.frame(myFilter) 
  myFilter = myFilter %>% mutate(geneId = colnames(myFilter)) %>% 
    pivot_longer(-geneId) %>% filter(value > 1)
  identicalARG = identicalARG %>% 
    filter(gene1 %in% myFilter$geneId, gene2 %in% myFilter$geneId)
  
  #Build the graph again and evaluate cliques
  identicalARG = graph_from_data_frame(data.frame(
    from = identicalARG$gene1, to = identicalARG$gene2
  ), directed = F)
  
  
  identicalARG = map_df(sapply(max_cliques(identicalARG), names), function(x){
    data.frame(geneId = x)
  }, .id = "ARGgroup")
  
  for(i in unique(identicalARG$ARGgroup)[n_distinct(identicalARG$ARGgroup):2]){
    identicalARG = identicalARG %>% 
      filter(ARGgroup == i | 
               (ARGgroup != i & !geneId %in% geneId[ARGgroup == i]))
  }
  
  #Per ARG group, filter out the best by using the gene statistics
  identicalARG = identicalARG %>% 
    left_join(genesDetected %>% 
                mutate(geneId = as.character(geneId)) %>%
                select(geneId, gene, subtype, startPerc, startDepth, cover1), 
              by = "geneId") %>% 
    mutate(val = startPerc * startDepth * cover1) %>% 
    group_by(ARGgroup) %>% mutate(keep = val == max(val)) %>% ungroup()
  
  #Only keep the genes that are the best in their cluster (i.e. remove duplicates)
  result = result %>% 
    left_join(identicalARG %>% select(geneId, keep), by = "geneId") %>% 
    mutate(keep = replace_na(keep, T)) %>% 
    filter(keep) %>% select(-keep)
  
  #---- Collapse overlapping bacteria ---
  #--------------------------------------
  bactGroups = blastOut %>% 
    rowwise() %>% 
    mutate(
      ident = identity / align_len,
      coverage = min(align_len / query_len, 1),
    ) %>%
    filter(coverage >= 0.9, ident >= 0.95) %>%
    group_by(segmentId, geneId, hitId, accession) %>% 
    filter(bit_score == max(bit_score)) %>% #Remove multiple hits (keep best)
    slice(1) %>% ungroup() %>% 
    select(segmentId, geneId, taxid, accession, plasmid, genus, species, strain, 
           extra, bit_score, coverage, ident) %>% 
    left_join(pathData %>% 
                select(segmentId, order, startOrientation, dist, depth, type) %>% 
                distinct(), by = "segmentId") %>% 
    filter(!is.na(dist), order < 3) %>% 
    left_join(genesDetected %>% 
                mutate(geneId = as.character(geneId)) %>% 
                select(geneId, gene, subtype), by = "geneId") %>% 
    rowwise() %>% 
    mutate(pathScore = bit_score * coverage / (max(dist,0) * 0.01 + 1)) %>% 
    group_by(geneId, genus, species) %>% 
    filter(n_distinct(startOrientation[order == 2]) == 2 | any(type == "fragment")) %>% 
    # group_by(segmentId) %>% filter(pathScore == max(pathScore)) %>%
    ungroup() 

  bactGroups = map_df(unique(bactGroups$geneId), function(x){
    x = bactGroups %>% filter(geneId == x) %>% 
      pull(taxid) %>% unique()
    if(length(x) > 1) {
      combn(x, 2) %>% t() %>% as.data.frame()
    } else {
      data.frame()
    }
    
  }) 
  
  #Similar to the ARG clustering, use graphs / cliques to detect overlap
  if(nrow(bactGroups) > 0){
    bactGroups = graph_from_data_frame(bactGroups  %>% distinct() %>% 
                                         select(from = V1, to = V2), directed = F)
    
    bactGroups = map_df(sapply(max_cliques(bactGroups), names), function(x){
      data.frame(taxid = x)
    }, .id = "bactGroup")
    
    for(i in unique(bactGroups$bactGroup)[n_distinct(bactGroups$bactGroup):2]){
      bactGroups = bactGroups %>% 
        filter(bactGroup == i | 
                 (bactGroup != i & !taxid %in% taxid[bactGroup == i]))
    }
    bactGroups = bactGroups %>% mutate(bactGroup = as.integer(bactGroup))
    
    test = result %>% 
      left_join(bactGroups %>% mutate(taxid = as.integer(taxid)), by = "taxid")
  } else {
    test = result %>% 
      group_by(genus, species) %>% 
      mutate(bactGroup = cur_group_id())

  }
  
  #Only keep bacteria that are really present
  test = test %>% 
    filter(taxid %in% allBact2$taxid) %>% 
    group_by(genus, species) %>% 
    mutate(bactGroup = ifelse(is.na(bactGroup), max(.$bactGroup, na.rm = T) + 
                                cur_group_id(), bactGroup)) %>% 
    group_by(gene, genus) %>% 
    filter(pathScore == max(pathScore)) %>% slice(1)  %>% 
    group_by(bactGroup) %>%
    mutate(bactGroup = paste(paste(genus, species) %>% unique(), collapse = ", "))
  
  nodes = data.frame(
    id = c(test$gene, test$bactGroup) %>% unique(),
    label = c(test$gene, test$bactGroup) %>% unique()
  )
  
  edges = 
    data.frame(
      from = test$bactGroup,
      to = test$gene
    ) %>% distinct()
  
  visNetwork(nodes, edges, height = "1000px")
  
# })

# results = results %>% select(pipelineId, everything())
# write_csv(results, 'sideStuff/annotationResults_4.csv')

