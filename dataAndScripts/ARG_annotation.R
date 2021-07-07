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

softmax = function(vals, normalise = F, log = T){
  if(normalise) vals = vals / max(vals)
  if(log) vals = log(vals)
  return(exp(vals) / sum(exp(vals)))
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
i = 17
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
    #Extract subspecies, strain and plasmid info
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
    #Get the coverage and identity of each alignment
    rowwise() %>% 
    mutate(
      ident = identity / align_len,
      coverage = min(align_len / query_len, 1),
    ) %>% ungroup() 
  
  #Fix some known issues with naming
  blastOut[blastOut$genus == "Enterobacter" & 
             blastOut$species == "aerogenes","genus"] = "Klebsiella"
  
  
  #Get the KC and LN for the blasted segments
  segmentInfo =  read_csv(paste0(sample, "/blastSegments.csv"), 
                          col_types = cols()) %>%
    select(-sequence, -name, -geneId)
 
  #Get information about the blasted seq where they are in the path of GFA
  pathData = list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T)
  
  pathData = map_df(pathData, function(myGFA){
    #Read the GFA
    gfa = gfa_read(myGFA)
    
    #If there are linked segments ...
    if(nrow(gfa$links) > 0){
      
      #Get the (largest) start segment
      segmentOfInterest = gfa$segments %>% 
        filter(str_detect(name, "_start$")) %>%
        filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% 
        slice(1) %>% pull(name)
      
      #Get all semenents in path to start
      x = gfa_pathsToSegment(gfa, segmentOfInterest, returnList = T, pathSegmentsOnly = T) %>% 
        map_df(function(path){
          data.frame(
            pathId = path$id,
            startOrientation = path$startOrientation,
            segmentId = path$segmentOrder)
        })
      
      #If no paths, return empty frame
      if(nrow(x) == 0 ){
        data.frame()
      }
      #Otherwise add the path data 
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
            #These files contain full GFA structures, no fragments
            type = "full")
      }
      
      
      
    } else if(nrow(gfa$segments) > 0){
      #Case where there are no links, just one or more segments
      gfa$segments %>% 
        select(segmentId = name, KC, LN) %>% 
        mutate(
          pathId = 1:n(),
          startOrientation = 0,
          geneId = str_extract(segmentId, "^\\d+"),
          dist = ifelse(str_detect(segmentId, "_start$"), -LN, 0),
          order = 1,
          type = "full"
        )
      
    } else{
      #File is empty
      data.frame()
    }
    
  })
  
  #Also add the fragment data to the paths
  fragments = gfa_read(paste0(sample, "/fragmentsOnly.gfa"))

  if(nrow(fragments$segments) > 0){
    #Fragments were merged, so treated as a non-start piece
    fragments = fragments$segments %>% 
      select(segmentId = name, LN, KC, geneId) %>% 
      group_by(geneId) %>% 
      mutate(pathId = 1, startOrientation = 0:(n()-1),
             order = 1, type = "fragment", dist = -1)
  } else {
    fragments = data.frame()
  }
  
  pathData = bind_rows(pathData, fragments)  %>%  ungroup()
  
  
  pathData = pathData %>% group_by(geneId, pathId, startOrientation) %>% 
    mutate(
      depth = KC / LN,
      totalLN = sum(LN[dist >= 0])) %>% ungroup()
  
  #---- Remove duplicate genes ---
  #-------------------------------
  
  #Get the segment info for ALL segments (not just blasted ones) 
  getSegments = map_df(list.files(
    paste0(sample, "/genesDetected/simplifiedGFA"), 
    pattern = ".gfa", full.names = T), function(x) gfa_read(x)$segments)
  
  #Only keep the start segments and the ones directly connecting to them 
  myFilter = pathData %>% filter(order < 3) %>% pull(segmentId) %>% unique()
  getSegments = getSegments %>% filter(name %in% myFilter)
  
  #Write these segments to a fasta file
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
  #Usearch outputs the pairs that are above identity threshold
  identicalARG = read_delim(sprintf("%s/ARGsimilarities.out", sample), "\t", col_names = F) %>% 
    filter(X1 != X2) %>% 
    mutate(X1 = str_remove(X1, "_RC$"), X2 = str_remove(X2, "_RC$")) %>% 
    distinct() %>% rowwise() %>% 
    #Create a group var of the gene pairs and add start info
    mutate(
      start = str_detect(X1, "_start") & str_detect(X2, "_start"),
      groupId = paste(sort(c(str_extract(X1, "^\\d+"),
                             str_extract(X2, "^\\d+"))), collapse = "_")) %>% 
    #Add the path data
    left_join(pathData %>% filter(order < 3) %>% 
                select(segmentId, d1 = depth, o1 = order, LN1 = LN, KC1 = KC) %>% 
                distinct(), by = c("X1" = "segmentId")) %>% 
    left_join(pathData %>% filter(order < 3) %>% 
                select(segmentId, d2 = depth, o2 = order, LN2 = LN, KC2 = KC) %>% 
                distinct(), by = c("X2" = "segmentId")) %>% 
    #Only keep groups that at least have overlapping start segment
    group_by(groupId) %>% 
    filter(any(start)) %>% ungroup() %>% 
    select(-start, -groupId) %>% distinct() %>%
    #Only keep segments that match in same distance (e.g. start vs start)
    filter(o1 == o2) %>% 
    rowwise() %>% 
    mutate(
      gene1 = str_extract(X1, "^\\d+"),
      gene2 = str_extract(X2, "^\\d+"),
      LNratio = min(LN1, LN2) / max(LN1, LN2),
      LNoverlap = min(LN1, LN2)
    ) %>% 
    #Remove alignments that are too big in size difference, unless start
    filter(LNratio > 0.9 | (X3 == 0 & o1 == 1 & LNoverlap > 250)) %>% 
    mutate(gene = ifelse(d1 > d2, gene1, gene2),
           depth = ifelse(d1 > d2, d1, d2)) %>% ungroup()
  
  #Build a graph to see how the genes are connected
  myFilter = graph_from_data_frame(data.frame(
    from = identicalARG$gene1, to = identicalARG$gene2
  ), directed = F)
  
  #Only consider similarities if at least two pieces align
  myFilter = as_adjacency_matrix(myFilter, sparse = F)
  myFilter = as.data.frame(myFilter) 
  myFilter = myFilter %>% mutate(geneId = colnames(myFilter)) %>% 
    pivot_longer(-geneId) %>% filter(value > 0)
  identicalARG = identicalARG %>% 
    filter(gene1 %in% myFilter$geneId, gene2 %in% myFilter$geneId)
  
  #Build the graph again and evaluate cliques
  identicalARG = graph_from_data_frame(data.frame(
    from = identicalARG$gene1, to = identicalARG$gene2
  ), directed = F)
  
  
  identicalARG = map_df(sapply(max_cliques(identicalARG), names), function(x){
    data.frame(geneId = x)
  }, .id = "ARGgroup")
  
  #Starting with the largest clique, make them unique by removing 
  #duplicates in other ones
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
    mutate(val = startPerc * startDepth * cover1,
           ARGgroup = as.integer(ARGgroup)) %>% 
    group_by(ARGgroup) %>% mutate(keep = val == max(val)) %>% ungroup()
  

  # ---- Filter / group bacteria ----
  #**********************************
  
  allBact = blastOut %>%
    select(segmentId, geneId, hitId, bit_score, coverage,
           accession, taxid, genus, species, extra, plasmid) %>%
    group_by(segmentId, geneId, hitId, accession) %>%
    filter(bit_score == max(bit_score)) %>% 
    slice(1) %>% ungroup() %>%
    #Add the path data
    left_join(pathData %>% 
                select(segmentId, order, startOrientation, dist) %>% 
                mutate(startOrientation = ifelse(str_detect(segmentId, "_start$"),
                                                 -1, startOrientation)) %>% 
                distinct(), by = "segmentId") %>% 
    #Remove those not in a path
    filter(!is.na(order)) %>% rowwise() %>% 
    mutate(pathScore = bit_score * coverage / (max(dist,0) * 0.01 + 1)) %>% 
    # filter(order < 4) %>%
    group_by(geneId, accession, taxid, order, startOrientation) %>% 
    #Get the best paths per accession
    filter(pathScore == max(pathScore)) %>% slice(1) %>% 
    group_by(geneId, accession, startOrientation) %>% 
    filter(all(2:max(order) %in% order) | all(order == 1)) %>%
    group_by(geneId, accession, taxid, genus, species, plasmid) %>% 
    summarise(n = n(), pathScore = sum(pathScore), 
              startOrientation = paste(unique(startOrientation), collapse = ","),
              orders = paste(order, collapse = ","), .groups = "drop") %>% 
    left_join(genesDetected %>%
                mutate(geneId = as.character(geneId)) %>%
                select(geneId, gene, subtype)) 

  #Filter bacteria based off the ARG groups
  bactGroups = allBact %>% 
    left_join(identicalARG %>% select(geneId, ARGgroup), by = "geneId") %>% 
    group_by(geneId) %>%
    mutate(ARGgroup = ifelse(is.na(ARGgroup), cur_group_id() + 
                               max(.$ARGgroup, na.rm = T), ARGgroup)) %>% 
    group_by(ARGgroup) %>% 
    filter(pathScore >= quantile(pathScore, 0.75)) %>% 
    ungroup()
  
  result = bactGroups %>% 
    select(ARGgroup, pathScore, taxid, genus, species) %>% 
    distinct()
  
  myFilter = result %>% select(taxid, genus, species, ARGgroup) %>% distinct() %>% 
    group_by(taxid) %>% mutate(n = n()) %>% 
    group_by(genus, species, ARGgroup) %>% filter(n == max(n))
  
  
  result = result %>% 
    filter(taxid %in% myFilter$taxid) %>% 
    group_by(ARGgroup, taxid) %>% 
    filter(pathScore == max(pathScore)) %>% ungroup() %>% 
    mutate(
      name = paste(str_extract(genus, "^.."), str_extract(species, "^.."), taxid)
    ) %>% 
    group_by(ARGgroup) %>%
    filter(pathScore >= quantile(pathScore, 0.9)) %>% 
    ungroup() %>% 
    select(from = name, to = ARGgroup, weight = pathScore)
  
  myGraph = graph_from_data_frame(result, directed = F)
  
  myGraph = data.frame(
    from = names(components(myGraph)$membership),
    membership = components(myGraph)$membership
  )
  
  result = result %>% left_join(myGraph, by = "from") %>% 
    group_by(membership, from) %>% 
    mutate(val = sum(weight)) %>% ungroup()
  
  bact = result %>% group_by(from, membership) %>% 
    summarise(val = max(val)) %>% 
    group_by(membership) %>% 
    mutate(prob = softmax(val, T, T)) %>% 
    arrange(membership, desc(prob))
  
  myGenes = bactGroups %>% select(geneId, ARGgroup) %>% 
    distinct() %>% left_join(
      genesDetected %>% mutate(geneId = as.character(geneId)), by = "geneId")
  
  
  genes = result %>% select(ARGgroup = to, membership) %>% 
    distinct() %>% left_join(
      myGenes %>% 
        mutate(
          val = startPerc * startDepth * cover1
        ) %>% 
        group_by(ARGgroup) %>% 
        filter(val == max(val)) %>% 
        select(ARGgroup, gene, subtype, cover1, type),
      by = "ARGgroup"
    ) %>% arrange(membership, gene)
  
  #PLOT
  if(F){
    test = result %>% 
      select(species, ARGgroup, bactGroup) %>% 
      distinct()
    
    nodes = data.frame(
      id = c(unique(test$species), unique(test$ARGgroup)),
      label = c(unique(test$species), unique(test$ARGgroup)),
      shape = c(rep("square", n_distinct(test$species)),
                rep("circle", n_distinct(test$ARGgroup)))
    )
    
    edges = data.frame(
      from = test$species,
      to = test$ARGgroup
    )
    
    visNetwork(nodes, edges, height = "1200px") %>% 
      visPhysics(stabilization = list(iterations = 5),
                 solver = "repulsion")
    
  }
 
  sampleName
# })

  