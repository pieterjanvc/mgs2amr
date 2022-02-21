#**************************
# ---- ARG ANNOTATION ----
#*************************

suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(visNetwork))
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

args = commandArgs(trailingOnly = TRUE) 
baseFolder = formatPath(args[[1]], endWithSlash = T)
database = args[[2]]
runId = as.integer(args[[3]])
verbose = abs(as.integer(args[[4]]))
pipelineIds = str_trim(unlist(strsplit(args[[5]], ",")))
generateReport = as.logical(args[[6]])
# forceRedo = as.logical(as.logical(args[[7]]))
forceRedo = T
  
maxCPU = 14

#Load the ARG and the sample list to process
myConn = myConn = dbConnect(SQLite(), database,  synchronous = NULL)

ARG = dbReadTable(myConn, "ARG") 

#Pick the samples to process
if(forceRedo & length(pipelineIds) > 0){
  #Redo (finished) samples
  toProcess = dbReadTable(myConn, "pipeline") %>% 
    filter(pipelineId %in% local(pipelineIds))
  
} else if(length(pipelineIds) > 0) {
  #Only process specific unprocessed samples
  toProcess = dbReadTable(myConn, "pipeline") %>% 
    filter(statusCode == 4 & pipelineId %in% local(pipelineIds))
} else {
  #All unprocessed
  toProcess = dbReadTable(myConn, "pipeline") %>% 
    filter(statusCode == 4)
}


#Generate a list out of the settings file
settings = readLines(paste0(baseFolder, "settings.txt"))
settings = settings[str_detect(settings,"^\\s*[^=#]+=.*$")]
settings = str_match(settings, "\\s*([^=\\s]+)\\s*=\\s*(.*)")
settings = setNames(str_trim(settings[,3]), settings[,2])

#Get the AMR prediction models
predictionModels = readRDS(
  sprintf("%sdataAndScripts/predictionModels.rds",baseFolder))
myAntibiotics = tbl(myConn, "antibiotics") %>% 
  filter(name %in% local(names(predictionModels))) %>% 
  as.data.frame()

bactGenomeSize = 
  read.csv(sprintf("%sdataAndScripts/bactGenomeSize.csv",baseFolder))

dbDisconnect(myConn)

#Check if the temp folders still exist
toProcess$exists = dir.exists(toProcess$tempFolder)
if(!all(toProcess$exists)){
  cat("\n!! The data for the following pipelineIds cannot be found", 
      "(deleted? moved?) and will be skipped:", 
      paste(toProcess$pipelineId[!toProcess$exists], collapse = ","), "\n\n")
  
  toProcess = toProcess %>% filter(exists)
}

# ---- FUNCTIONS ----
#********************
softmax = function(vals, normalise = F, log = T){
  if(normalise) vals = vals / max(vals)
  if(log) vals = log(vals)
  return(exp(vals) / sum(exp(vals)))
}

normalise = function(x){
  (x - min(x)) / (max(x) - min(x))
}

# ---- Annotate for all pipelineIds ----
#***************************************

options(readr.num_columns = 0)

#Check if there are any samples to process
if(nrow(toProcess) == 0) {
  
  newLogs = data.frame(
    timeStamp = as.integer(Sys.time()), 
    actionId = 15, 
    actionName = "Nothing to annotate",
    runId = runId,
    tool = "ARG_annotation.R"
    ) %>% select(runId,tool,timeStamp,actionId,actionName)
  
  myConn = myConn = dbConnect(SQLite(), database, synchronous = NULL)
  
  q = dbSendStatement(
    myConn, 
    "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
    params = unname(as.list(newLogs)))
  q = dbClearResult(q)
  
  dbDisconnect(myConn)
  
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                      "No new pipelines found to annotate. Exit script\n\n")}
  
} else {
  

  sampleIndex = 1:nrow(toProcess)
  registerDoParallel(cores=min(maxCPU, length(sampleIndex)))
  
  result = foreach(sampleIndex = sampleIndex) %dopar% {
    
    #Process each sample
    tryCatch({
      success = F
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
      
      # ---- FILTERING DATA ----
      #*************************
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Loading BLASTn output ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                    "Start loading BLASTn output"))
      
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
        filter(!genus %in% c("uncultured", "Uncultured", "mixed", "Bacterium"),
               !species %in% c("bacterium", "Bacterium") & 
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
      
      # blastOut[blastOut$genus == "Enterobacter","species"] = "sp."
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                    "Finished loading BLASTn output"))
      
      
      # #test
      # toRemove = ARG %>% filter(str_detect(gene, "tet.+\\/")) %>% 
      #   pull(geneId)
      # genesDetected = genesDetected %>% filter(!geneId %in% toRemove)
      # blastOut = blastOut %>% filter(!geneId %in% toRemove)
      
      # ---- MAPPING DATA TO GFA ----
      #******************************
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Mapping alignments to GFA ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                                    "Start mapping"))
      
      
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
            dplyr::slice(1) %>% pull(name)
          
          if(length(segmentOfInterest) == 0){
            return(data.frame())
          }
          
          #Get all semenents in path to start
          x = gfa_pathsToSegment(gfa, segmentOfInterest, returnList = T, 
                                 pathSegmentsOnly = T, verbose = F) %>% 
            map_df(function(path){
              data.frame(
                pathId = path$id,
                orientation = path$orientation,
                segmentId = path$segmentOrder)
            })
          
         
          #If no paths, return empty frame
          if(nrow(x) == 0 ){
            data.frame()
          }
          #Otherwise add the path data 
          else{
             
            x = x %>% 
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
            
            if(any(x$pathType == 2)){
              
              #In case of circular loop, only keep 1 path that has least competition
              toRemove = x %>% filter(pathType != 2) %>% group_by(pathId, orientation) %>% 
                summarise(LN = sum(LN) - 30*(n() - 1), .groups = "drop") %>% 
                filter(LN == max(LN)) %>% slice(1) %>% pull(orientation)
              
              x = x %>% filter(!(pathType == 2 & orientation == toRemove))
            }
            
            return(x)
          }
          
          
          
        } else if(nrow(gfa$segments) > 0){
          #Case where there are no links, just one or more segments
          gfa$segments %>% 
            select(segmentId = name, KC, LN) %>% 
            mutate(
              pathId = 1:n(),
              orientation = 0,
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
          mutate(pathId = 1, orientation = 0:(n()-1),
                 order = 1, type = "fragment", dist = -1)
      } else {
        fragments = data.frame()
      }
      
      pathData = bind_rows(pathData, fragments)  %>%  ungroup()
      
      
      pathData = pathData %>% group_by(geneId, pathId, orientation) %>% 
        mutate(
          depth = KC / LN,
          totalLN = sum(LN[dist >= 0])) %>% ungroup()
      
      #Only keeps paths that have results from blast
      pathData = pathData %>% left_join(data.frame(
        segmentId = unique(blastOut$segmentId),
        result = T
      ),by = "segmentId") %>% 
        group_by(geneId, pathId) %>% 
        filter(!any(is.na(result))) %>% ungroup() %>% select(-result)
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, 
                                    "Finished mapping"))
      
      
      #---- Remove duplicate genes ---
      #-------------------------------
      
      # if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
      #                     "Detect duplicate genes ... ")}
      # newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, 
      #                               "Start detecting duplicate genes"))
      # 
      # #Get the segment info for ALL segments (not just blasted ones) 
      # getSegments = map_df(list.files(
      #   paste0(sample, "/genesDetected/simplifiedGFA"), 
      #   pattern = ".gfa", full.names = T), function(x) gfa_read(x)$segments)
      # 
      # #Only keep the start segments and the ones directly connecting to them 
      # myFilter = pathData %>% filter(order < 3) %>% pull(segmentId) %>% unique()
      # getSegments = getSegments %>% filter(name %in% myFilter)
      # 
      # #Write these segments to a fasta file
      # fasta_write(getSegments$sequence, 
      #             sprintf("%s/ARGsim.fasta", sample),
      #             getSegments$name, type = "n")
      # 
      # #Add the reverse complement (usearch does not do that when calc_distmx)
      # system(sprintf("%s -fastx_revcomp %s -label_suffix _RC -fastaout %s >/dev/null 2>&1",
      #                settings["usearch"],
      #                sprintf("%s/ARGsim.fasta", sample),
      #                sprintf("%s/ARGsim_RC.fasta", sample)))
      # system(sprintf("cat %1$s %2$s > %3$s; rm %1$s %2$s >/dev/null 2>&1",
      #                sprintf("%s/ARGsim.fasta", sample),
      #                sprintf("%s/ARGsim_RC.fasta", sample),
      #                sprintf("%s/ARGsimilarities.fasta", sample)))
      # 
      # #Use cluster_fast to reduce number of segments by grouping in identity clusters
      # system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.1 %s",
      #                settings["usearch"],
      #                sprintf("%s/ARGsimilarities.fasta", sample),
      #                sprintf("%s/ARGsimilarities.out", sample),
      #                ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))
      # 
      # #Group identical ARG and remove the duplicate (weak) ones
      # #Usearch outputs the pairs that are above identity threshold
      # identicalARG = read.delim(sprintf("%s/ARGsimilarities.out", sample), 
      #                           header = F, sep = "\t") %>% 
      #   rename(X1 = V1, X2 = V2, X3 = V3) %>% 
      #   filter(X1 != X2) %>% 
      #   mutate(X1 = str_remove(X1, "_RC$"), X2 = str_remove(X2, "_RC$")) %>% 
      #   distinct() %>% rowwise() %>% 
      #   #Create a group var of the gene pairs and add start info
      #   mutate(
      #     start = str_detect(X1, "_start") & str_detect(X2, "_start"),
      #     groupId = paste(sort(c(str_extract(X1, "^\\d+"),
      #                            str_extract(X2, "^\\d+"))), collapse = "_")) %>% 
      #   #Add the path data
      #   left_join(pathData %>% filter(order < 3) %>% 
      #               select(segmentId, d1 = depth, o1 = order, LN1 = LN, KC1 = KC) %>% 
      #               distinct(), by = c("X1" = "segmentId")) %>% 
      #   left_join(pathData %>% filter(order < 3) %>% 
      #               select(segmentId, d2 = depth, o2 = order, LN2 = LN, KC2 = KC) %>% 
      #               distinct(), by = c("X2" = "segmentId")) %>% 
      #   #Only keep groups that at least have overlapping start segment
      #   group_by(groupId) %>% 
      #   filter(any(start)) %>% ungroup() %>% 
      #   select(-start, -groupId) %>% distinct() %>%
      #   #Only keep segments that match in same distance (e.g. start vs start)
      #   filter(o1 == o2) 
      # 
      # if(nrow(identicalARG) > 0){
      #   identicalARG = identicalARG %>% 
      #     rowwise() %>% 
      #     mutate(
      #       gene1 = str_extract(X1, "^\\d+"),
      #       gene2 = str_extract(X2, "^\\d+"),
      #       LNratio = min(LN1, LN2) / max(LN1, LN2),
      #       LNoverlap = min(LN1, LN2)
      #     ) %>% 
      #     #Remove alignments that are too big in size difference, unless start
      #     filter(LNratio > 0.9 | (X3 == 0 & o1 == 1 & LNoverlap > 250)) %>% 
      #     mutate(gene = ifelse(d1 > d2, gene1, gene2),
      #            depth = ifelse(d1 > d2, d1, d2)) %>% ungroup()
      #   
      #   #Build a graph to see how the genes are connected
      #   myFilter = graph_from_data_frame(data.frame(
      #     from = identicalARG$gene1, to = identicalARG$gene2
      #   ), directed = F)
      #   
      #   #Only consider similarities if at least two pieces align
      #   myFilter = as_adjacency_matrix(myFilter, sparse = F)
      #   myFilter = as.data.frame(myFilter) 
      #   myFilter = myFilter %>% mutate(geneId = colnames(myFilter)) %>% 
      #     pivot_longer(-geneId) %>% filter(value > 0)
      #   identicalARG = identicalARG %>% 
      #     filter(gene1 %in% myFilter$geneId, gene2 %in% myFilter$geneId)
      #   
      #   #Build the graph again and evaluate cliques
      #   identicalARG = graph_from_data_frame(data.frame(
      #     from = identicalARG$gene1, to = identicalARG$gene2
      #   ), directed = F)
      #   
      #   
      #   identicalARG = map_df(sapply(max_cliques(identicalARG), names), function(x){
      #     data.frame(geneId = x)
      #   }, .id = "ARGgroup")
      #   
      #   #Starting with the largest clique, make them unique by removing 
      #   #duplicates in other ones
      #   for(i in unique(identicalARG$ARGgroup)[n_distinct(identicalARG$ARGgroup):2]){
      #     identicalARG = identicalARG %>% 
      #       filter(ARGgroup == i | 
      #                (ARGgroup != i & !geneId %in% geneId[ARGgroup == i]))
      #   }
      #   
      #   #Per ARG group, filter out the best by using the gene statistics
      #   identicalARG = identicalARG %>% 
      #     left_join(genesDetected %>% 
      #                 mutate(geneId = as.character(geneId)) %>%
      #                 select(geneId, startPerc, startDepth, cover), 
      #               by = "geneId") %>% 
      #     mutate(val = startPerc * startDepth * cover,
      #            ARGgroup = as.integer(ARGgroup)) %>% 
      #     group_by(ARGgroup) %>% mutate(keep = val == max(val)) %>% ungroup()
      # } else {
      #   identicalARG = data.frame(geneId = as.character(), 
      #                             ARGgroup = as.numeric(), keep = as.logical())
      # }
      # 
      # 
      # 
      # if(verbose > 0){cat("done\n")}
      # newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, 
      #                               "Finished detecting duplicate genes"))
      # 
      # 
      # #Add the grouping info to detected genes
      # genesDetected = genesDetected %>%
      #   left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>%
      #   select(pipelineId, geneId, gene, subtype, everything())
      # 
      # genesDetected = genesDetected %>% 
      #   select(-matches("ARGgroup")) %>% 
      #   left_join(identicalARG %>% 
      #               select(geneId, ARGgroup, keep) %>% 
      #               mutate(geneId = as.integer(geneId)), by = "geneId") %>% 
      #   mutate(keep = ifelse(is.na(keep), T, keep)) 
      # 
      # genesDetected[is.na(genesDetected$ARGgroup),] = 
      #   genesDetected[is.na(genesDetected$ARGgroup),] %>% 
      #   group_by(geneId) %>% 
      #   mutate(ARGgroup = cur_group_id() + 
      #            max(c(0, genesDetected$ARGgroup), na.rm = T)) %>% 
      #   ungroup() 
      
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
      
      blastSegments = 
        read.csv(paste0(sample, "/blastSegments.csv")) %>% 
        select(name, LN, KC)
      blastOut = blastOut %>% left_join(blastSegments, 
                                        by = c("segmentId" = "name")) %>% 
        mutate(KC = KC * coverage * ident)
      
      # blastOut = blastOut %>% group_by(clusterId) %>% 
      #   filter(bit_score >= 0.90 * max(bit_score)) %>% ungroup()
      
      #Make sure that taxid always has the same genus / species name and vice versa
      blastOut = blastOut %>% 
        # mutate(bact = paste(genus, species)) %>% 
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
      
      allBact = blastOut %>%
        # filter(geneId == "3323") %>%
        select(segmentId, geneId, bit_score, coverage,
               accession, taxid, genus, species, extra, plasmid, KC, LN) %>%
        group_by(segmentId, geneId, accession) %>%
        filter(bit_score == max(bit_score)) %>% 
        dplyr::slice(1) %>% ungroup() %>%
        mutate(start = str_detect(segmentId, '_start$')) %>%
        
        mutate(pathScore = bit_score * coverage) %>% 
        group_by(segmentId, taxid, plasmid) %>% 
        filter(pathScore == max(pathScore)) %>% dplyr::slice(1) %>% ungroup() %>% 
        left_join(pathData %>% 
                    select(pathId, segmentId, order, orientation) %>% 
                    mutate(orientation = ifelse(str_detect(segmentId, "_start$"),
                                                     -1, orientation)) %>% 
                    distinct(), by = "segmentId") %>% 
        filter(!is.na(order))  %>% 
        group_by(geneId, pathId, taxid, plasmid, orientation) %>% 
        filter(all(2:max(order) %in% order) | all(order == 1)) %>% 
        group_by(geneId, taxid, plasmid, orientation, pathId) %>% 
        mutate(fullPath = sum(pathScore)) %>% 
        group_by(geneId, taxid, plasmid, orientation) %>% 
        filter(pathId == pathId[fullPath == max(fullPath)][1]) %>% 
        group_by(geneId, taxid, plasmid) %>% 
        mutate(fullPath = sum(pathScore)) %>% 
        group_by(geneId, taxid) %>% 
        filter(plasmid == plasmid[fullPath == max(fullPath)][1])
      
      
      #Remove genes with no blast data
      geneIds = genesDetected %>% pull(geneId)
      geneIds = geneIds[geneIds %in% unique(allBact$geneId)]
      
      #Check the overlap between segments of different species
      # if one is completely contained within another (and smaller), it's a duplicate
      allBact = map_df(geneIds, function(geneId){
        
        #Examine the paths of one ARG
        myGene = allBact %>% filter(geneId %in% {{geneId}})
        
        #Get all the segmetents per path
        bactSegments = myGene %>% 
          group_by(segmentId, taxid, pathId) %>% 
          summarise(.groups = "drop") %>% 
          group_by(taxid, pathId) %>% 
          summarise(x = list(segmentId), .groups = "drop") 
        
        if(length(bactSegments$taxid) == 0){
          
          #There is no blast output for this gene
          return(data.frame())
          
        } else if(length(bactSegments$taxid) == 1){
          
          #There is only one taxid for the gene
          myGene = myGene %>% 
            group_by(geneId, taxid, genus, species, plasmid) %>% 
            summarise(n = n(), 
                      extension = sum(pathScore[!start]),
                      pathScore = sum(pathScore), 
                      KC = sum(KC), LN = sum(LN),
                      .groups = "drop") %>%
            group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
            arrange(desc(pathScore))
          
          myTaxid = bactSegments$taxid
          
          return(myGene %>% filter(taxid %in% myTaxid) %>% 
                   group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
                   dplyr::slice(1) %>% ungroup())
          
        } else {
          
          #There are multiple taxids for the gene
          #Check per path if there are redundant taxids (contained within other)
          pathResult = map_df(unique(bactSegments$pathId), function(pathId){
            
            #Get segments for the path
            pathSegments = bactSegments %>% filter(pathId == {{pathId}})
            
            if(nrow(pathSegments) > 1){
              #Build a matrix for species overlap per segment
              gfaMatrix = sapply(pathSegments$x, function(x){
                len = length(x)
                sapply(pathSegments$x, function(y){
                  sum(x %in% y) / len
                })
              })
              
              rownames(gfaMatrix) = pathSegments$taxid
              colnames(gfaMatrix) = pathSegments$taxid
            }
            
            myPath = myGene %>% 
              filter(taxid %in% pathSegments$taxid) %>% 
              group_by(geneId, taxid, genus, species, plasmid) %>% 
              summarise(n = n(), 
                        extension = sum(pathScore[!start]),
                        pathScore = sum(pathScore), 
                        KC = sum(KC), LN = sum(LN),
                        .groups = "drop") %>%
              group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
              arrange(desc(pathScore))
            
            myTaxid = myPath$taxid %>% unique()
            i = 1
            
            while(i < length(myTaxid)){
              myIndex = which(colnames(gfaMatrix) == myTaxid[i])
              duplicates = data.frame(
                name = colnames(gfaMatrix),
                row = gfaMatrix[myIndex,],
                col = gfaMatrix[,myIndex]
              ) %>% filter(row == 1, col != 1) %>% pull(name)
              
              myTaxid = setdiff(myTaxid, duplicates)
              i = i + 1
            }
            
            return(myPath %>% filter(taxid %in% myTaxid) %>% 
                     group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
                     dplyr::slice(1) %>% ungroup())
            
          })
          
          return(pathResult)
          
        } 
        
      }) %>% mutate(depth = KC / LN) %>% 
        left_join(
          genesDetected %>% mutate(geneId = as.character(geneId)) %>% 
            select(geneId, gene, subtype, cover, type), 
          by = "geneId"
        )
      
      #Adjust the scored for genes by presences of other genes in the same bact
      # myData = allBact %>% filter(!geneId %in% AMRclusters$geneId)
      # bactGroupStart = max(AMRclusters$bactGroup)
      
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
            group_by(geneId, bact = taxid) %>%
            summarise(pathScore = round(max(pathScore),2), 
                      extension = round(max(extension),2), 
                      .groups = "drop") %>%
            left_join(genesDetected %>%
                        select(geneId, ARGgroup) %>%
                        mutate(geneId = as.character(geneId)), by = "geneId")
          
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
            geneId = colnames(geneMatrix)
          ) %>% left_join(
            genesDetected %>% select(geneId, cover) %>%
              mutate(geneId = as.character(geneId)), by = "geneId")
          
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
          myClusters = myClusters %>% group_by(bactGroup) %>% filter(!all(!primary)) %>% 
            ungroup() %>% mutate(bactGroup = bactGroupStart + as.factor(bactGroup) %>% as.integer()) %>% 
            pivot_longer(c(-taxid, -bactGroup, -primary), names_to = "geneId", values_to = "val") %>% 
            mutate(geneId = as.integer(geneId), taxid = as.integer(taxid)) #%>% 
          
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
            myData %>% group_by(taxid) %>% 
              summarise(
                depth = sum(depth * softmax(cover, normalise = T)),
                pathScore = sum(pathScore), 
                .groups = "drop"
              ),
            by = "taxid"
          ) 
        
        return(list(myClusters = myClusters, bactList = bactList))
      }
      
      #Get all the ARG that are exclusively found in genomes
      genomeARG = allBact %>% 
        group_by(geneId) %>% filter(all(!plasmid)) %>% 
        ungroup()
      
      #Check if there are any, proceed accordingly
      if(nrow(genomeARG) > 0){
        
        #Adjust the bact presence across ARG
        genomeARG = adjustBact(genomeARG)
        
        #Get the clusters and bact list
        AMRclusters = genomeARG$myClusters %>%
          select(-primary) %>% filter(val > 0) %>% 
          mutate(origin = "genome")
        
        bactList = genomeARG$bactList
        
        #Check if any of the remaining ARG match a genome ARG species
        genomeWithPlasmid = allBact %>% 
          filter(!geneId %in% genomeARG$myClusters$geneId) %>% 
          left_join(
            genomeARG$bactList %>% 
              # filter(val == 1) %>% #maybe remove?
              select(taxid, bactGroup, prob, genomeDepth = depth,
                     genomePathscore = pathScore, val) %>%  
              mutate(genomePathscore = genomePathscore),
            by = "taxid"
          ) %>% left_join(
            AMRclusters %>% 
              # filter(val == 1) %>% 
              select(taxid, geneId) %>% 
              left_join(genesDetected %>% 
                          select(geneId, cover, startDepth, type),
                        by = "geneId") %>% 
              group_by(taxid) %>% 
              summarise(
                genomeCover = weighted.mean(cover, cover),
                genomeType = case_when(
                  all(type == "noFragments") ~ "noFragments",
                  all(type == "fragmentsOnly") ~ "fragmentsOnly",
                  TRUE ~ "mixed"
                ), .groups = "drop"),
            by = "taxid"
          )
        
        #Only match if the quality of ARG is comparable (type, cover, ...)
        # EXPERIMENTAL ...
        genomeWithPlasmid = genomeWithPlasmid %>% 
          #Pathscore should be at in the top range
          group_by(geneId) %>% filter(pathScore >= 0.75*max(pathScore)) %>% 
          ungroup() %>% 
          filter(cover >= genomeCover*0.75, depth >= 0.25*genomeDepth,
                 genomeCover >= cover*0.75)
        
        #Only continue if there are matches
        if(nrow(genomeWithPlasmid) > 0){
          
          genomeWithPlasmid = bind_rows(
            genomeWithPlasmid %>% filter(is.na(bactGroup)),
            genomeWithPlasmid %>% 
              filter(!is.na(bactGroup)) %>% 
              group_by(geneId, bactGroup) %>%
              filter(prob == max(prob)) %>% dplyr::slice(1) %>% 
              ungroup()
          ) %>% 
            mutate(across(c(genomeDepth, genomePathscore, val), function(x) ifelse(is.na(x), 0, x)))
          
          #Keep the best scoring bacterium per gene
          genomeWithPlasmid = genomeWithPlasmid %>% 
            group_by(geneId) %>% filter(pathScore == max(pathScore))
          
          #Add the ARG to the exisintg clusters (no new clusters)
          AMRclusters = bind_rows(
            AMRclusters,
            genomeWithPlasmid %>% 
              filter(!is.na(bactGroup)) %>% 
              select(taxid, bactGroup, geneId) %>% 
              mutate(geneId = as.integer(geneId), origin = "linkedPlasmid")
          )
          
        }
        
        #Check if if there are any remaining ARG that are isolated plasmids
        isolatePlasmid = adjustBact(
          allBact %>% 
            filter(!geneId %in% AMRclusters$geneId), 
          bactGroupStart = max(AMRclusters$bactGroup))
        
        #Update the cluster info if any results
        if(nrow(isolatePlasmid$myClusters) > 0){
          AMRclusters = bind_rows(
            AMRclusters,
            isolatePlasmid$myClusters %>%
              select(-primary) %>% filter(val > 0) %>% 
              mutate(origin = "isolatePlasmid")
          )
          
          bactList = bind_rows(bactList, isolatePlasmid$bactList)
        }
        
      } else {
        
        #No pure genome ARG are found, only plasmids
        isolatePlasmid = adjustBact(allBact)
        
        AMRclusters = isolatePlasmid$myClusters %>%
          select(-primary) %>% filter(val > 0) %>% 
          mutate(origin = "isolatePlasmid")
        
        bactList = isolatePlasmid$bactList
      }

      #add the most likely gene match to the detectedARG
      genesDetected = genesDetected %>% select(-bactGroup) %>% left_join(
        AMRclusters %>% filter(val == 1 | is.na(val)) %>% 
          select(geneId, bactGroup) %>% distinct(),
        by = "geneId"
      )
      
      #Add the IDs
      AMRclusters = AMRclusters %>% mutate(
        pipelineId = {{myPipelineId}}, runId = {{runId}}
      )
      
      
      #Estimate the relative abundance of the bacteria
      bactList = bactList %>% 
        mutate(estimatedAbundance = depth * size * 1e6 / inputfileBP,
               across(c(prob, estimatedAbundance), round, digits = 4),
               runId = {{runId}})
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, 
                                    "Finished grouping bacteria"))
      
      #need to edit the database to add new amrCluster table and check genedDetected
      
      
      # ---- MAKE PREDICTIONS ----
      #***************************
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Predicting AMR ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, 
                                    "Start AMR prediction"))
      
      myInputs = genesDetected %>%
        mutate(id = paste(pipelineId,  bactGroup, sep = "_")) %>%
        group_by(pipelineId, bactGroup, gene) %>%
        mutate(val = startPerc * startDepth * cover) %>% 
        filter(val == max(val)) %>% dplyr::slice(1) %>%
        ungroup() %>% 
        select(id, gene, cover) %>%
        pivot_wider(gene, names_from = id, values_from = cover)
      
      myInputs = map_df(sapply(predictionModels, function(x) x$feature_names), function(y){
        data.frame(
          gene = y
        )
      }, .id = "antibiotic") %>%
        left_join(myInputs, by = "gene")
      
      myInputs[is.na(myInputs)] = 0
      
      predictions = map_df(myAntibiotics$name, function(myAntibiotic){
        
        testInput = myInputs %>% filter(antibiotic == myAntibiotic) %>%
          select(-antibiotic, -gene) %>% as.matrix() %>% t()
        
        data.frame(
          id = rownames(testInput),
          antibiotic = myAntibiotic,
          value = round(predict(predictionModels[[myAntibiotic]], testInput),4)
        ) %>% mutate(resistance = ifelse(value < 0.5, "S", 'R')) %>% 
          separate(id, into = c("pipelineId", "bactGroup"), "_")
        
      }) %>% left_join(myAntibiotics %>% select(name, antibioticId), 
                       by = c("antibiotic" = "name")) %>% 
        mutate(runId = runId)
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, 
                                    "Finished AMR prediction"))
      
      
      # ---- Save results ----
      #***********************
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Save all results / reports ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 12, 
                                    "Save results"))
      
      myConn = myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      
      #Add the ARGgroup to the detectedARG table (remove?)
      q = dbSendStatement(myConn, 
                    "UPDATE detectedARG 
            SET runId = ?, ARGgroup = ?, bactGroup = ?, plasmid = ?
            WHERE pipelineId = ? AND geneId = ?", 
                    params = list(genesDetected$runId, genesDetected$ARGgroup, 
                                  genesDetected$bactGroup, genesDetected$plasmid,
                                  genesDetected$pipelineId, genesDetected$geneId))
      q = dbClearResult(q)
      #Add the detected bacteria to the detectedBact table
      q = dbSendStatement(
        myConn, 
        sprintf("DELETE FROM detectedBact WHERE pipelineId = %i", 
                bactList$pipelineId[1]))
      q = dbClearResult(q)
      q = dbSendStatement(myConn, 
                    "INSERT INTO detectedBact (pipelineId, runId, taxId, bactGroup, 
            genus, species, prob, relAbundance, value) 
            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", 
                    params = bactList %>% 
                      select(pipelineId, runId, taxid, bactGroup, genus, species, 
                             prob, estimatedAbundance, pathScore) %>% 
                      as.list() %>% unname())
      q = dbClearResult(q)
      #Add the clusters to the AMRclusters table
      q = dbSendStatement(
        myConn, 
        sprintf("DELETE FROM AMRclusters WHERE pipelineId = %i", 
                bactList$pipelineId[1]))
      q = dbClearResult(q)
      q = dbSendStatement(myConn, 
                    "INSERT INTO AMRclusters (pipelineId, runId, bactGroup, 
            taxId, geneId, val, origin) 
            VALUES(?, ?, ?, ?, ?, ?, ?)", 
                    params = AMRclusters %>% 
                      select(pipelineId, runId, bactGroup, taxid, 
                             geneId, val, origin) %>% 
                      as.list() %>% unname())
      q = dbClearResult(q)
      #Add the AMT predictions bacteria to the AMRprediction table
      q = dbSendStatement(
        myConn, 
        sprintf("DELETE FROM AMRprediction WHERE pipelineId = %i", bactList$pipelineId[1]))
      q = dbClearResult(q)
      q = dbSendStatement(myConn, 
                    "INSERT INTO AMRprediction (pipelineId, runId, 
            bactGroup, antibioticId, resistance, value) 
            VALUES(?, ?, ?, ?, ?, ?)", 
                    params = predictions %>% 
                      select(pipelineId, runId, bactGroup, antibioticId,
                             resistance, value) %>% 
                      as.list() %>% unname())
      q = dbClearResult(q)
      dbDisconnect(myConn)
      # write_csv(genes, sprintf("%s/genes.csv", sample, sampleName))
      # write_csv(bact, sprintf("%s/bacteria.csv", sample, sampleName))
      
      #Render the markdown report if requested
      if(generateReport){
        # render(sprintf("%s/dataAndScripts/report.rmd", baseFolder),
        #        params = list(pipelineId = myPipelineId, 
        #                      baseFolder = baseFolder),
        #        output_file = sprintf("%s/%s_report.html", sample, sampleName))
      }
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 13, 
                                    "Annotation finished"))
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          sprintf("The pipeline (id %s) finished successfully\n\n",
                                  myPipelineId))}
      
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14, 
                                    "Pipeline finished successfully"))
      
      success = T
    },
    finally = {
      
      #Write the logs to the database, irrespective of success
      newLogs$runId = runId
      newLogs$tool = "ARG_annotation.R"
      newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
      
      myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      
      q = dbSendStatement(
        myConn, 
        "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
        params = unname(as.list(newLogs)))
      q = dbClearResult(q)
      
      #If pipeline successfully completed, note that in DB
      if(14 %in% newLogs$actionId){
      q = dbSendStatement(
          myConn, 
          "UPDATE pipeline SET statusCode = 5, inputFileBP = ?,
           statusMessage = 'Finished pipeline', modifiedTimestamp = ?
           WHERE pipelineId = ?",
          params = list(inputfileBP, as.character(Sys.time()),myPipelineId))
      q = dbClearResult(q)
      }
      
      dbDisconnect(myConn)
      cat(myPipelineId, ifelse(success, 'success', 'fail'), "\n")
    })
    
  }
}

