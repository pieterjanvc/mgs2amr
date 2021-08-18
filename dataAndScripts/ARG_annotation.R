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

library(parallel)

args = commandArgs(trailingOnly = TRUE) 
baseFolder = formatPath(args[[1]], endWithSlash = T)
runId = as.integer(args[[2]])
verbose = abs(as.integer(args[[3]]))
pipelineIds = str_trim(unlist(strsplit(args[[4]], ",")))
generateReport = as.logical(args[[5]])

#Get the arguments
# baseFolder = "/mnt/meta2amrData/meta2amr/"
# runId = 0
# verbose = 1
# pipelineIds = NULL
# generateReport = T

#Load the ARG and the sample list to process
myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))

ARG = dbReadTable(myConn, "ARG") 

toProcess = dbReadTable(myConn, "pipeline") %>% 
  filter(statusCode == 4)
if(length(pipelineIds) > 0){
  toProcess = toProcess %>% filter(pipelineId %in% pipelineIds)
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

bactGenomeSize = read_csv(sprintf("%sdataAndScripts/bactGenomeSize.csv",baseFolder))

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
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  q = dbExecute(
    myConn, 
    "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
    params = unname(as.list(newLogs)))
  dbDisconnect(myConn)
  
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                      "No new pipelines found to annotate. Exit script\n\n")}
  
} else {
  
  # UNCOMMENT if running in parallel
  
  # cl <- parallel::makeCluster(detectCores(), outfile = "")
  # x = clusterEvalQ(cl, {
  #   library(tidyverse)
  #   library(RSQLite)
  #   library(xgboost)
  #   library(igraph)
  #   library(gfaTools)
  # })
  # clusterExport(
  #   cl,
  #   varlist = c("baseFolder", "toProcess", "softmax", "settings", "ARG",
  #               "predictionModels", "myAntibiotics", "bactGenomeSize", "runId",
  #               "generateReport", "verbose"))
  # 
  # 
  # result = parLapply(cl, 1:nrow(toProcess), function(sampleIndex){
  
  for(sampleIndex in 1:nrow(toProcess)){
    
    #Process each sample
    tryCatch({
      
      sample = toProcess$tempFolder[sampleIndex]
      myPipelineId = toProcess$pipelineId[sampleIndex]
      sampleName = str_extract(sample, "[^\\/]+(?=_\\d+$)")
      
      inputfileBP = grep("sequences length",
                         readLines(sprintf("%s/metacherchant_logs/log", sample)),
                         value = T, fixed = T) %>%
        str_extract("([\\d'])+(?=\\s\\()") %>%
        str_replace_all("'", "") %>% as.numeric()
      
      #Grab the detected ARG from the previous step
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
      genesDetected = tbl(myConn, "detectedARG") %>% 
        filter(pipelineId == myPipelineId) %>% as.data.frame()
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
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                    "Finished loading BLASTn output"))
      
      
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
          
          #Get all semenents in path to start
          x = gfa_pathsToSegment(gfa, segmentOfInterest, returnList = T, 
                                 pathSegmentsOnly = T, verbose = F) %>% 
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
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, 
                                    "Finished mapping"))
      
      
      #---- Remove duplicate genes ---
      #-------------------------------
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Detect duplicate genes ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, 
                                    "Start detecting duplicate genes"))
      
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
      system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.1 %s",
                     settings["usearch"],
                     sprintf("%s/ARGsimilarities.fasta", sample),
                     sprintf("%s/ARGsimilarities.out", sample),
                     ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))
      
      
      #Group identical ARG and remove the duplicate (weak) ones
      #Usearch outputs the pairs that are above identity threshold
      identicalARG = read_delim(sprintf("%s/ARGsimilarities.out", sample), "\t", 
                                col_names = F, show_col_types = F) %>% 
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
                    select(geneId, startPerc, startDepth, cover1), 
                  by = "geneId") %>% 
        mutate(val = startPerc * startDepth * cover1,
               ARGgroup = as.integer(ARGgroup)) %>% 
        group_by(ARGgroup) %>% mutate(keep = val == max(val)) %>% ungroup()
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, 
                                    "Finished detecting duplicate genes"))
      
      
      #Add the grouping info to detected genes
      genesDetected = genesDetected %>% 
        left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>% 
        select(pipelineId, geneId, gene, subtype, everything())
      
      genesDetected = genesDetected %>% 
        select(-matches("ARGgroup")) %>% 
        left_join(identicalARG %>% 
                    select(geneId, ARGgroup, keep) %>% 
                    mutate(geneId = as.integer(geneId)), by = "geneId") %>% 
        mutate(keep = ifelse(is.na(keep), T, keep))
      
      genesDetected[is.na(genesDetected$ARGgroup),] = 
        genesDetected[is.na(genesDetected$ARGgroup),] %>% 
        group_by(geneId) %>% 
        mutate(ARGgroup = cur_group_id() + 
                 max(genesDetected$ARGgroup, na.rm = T)) %>% 
        ungroup()
      
      
      # ---- Filter / group bacteria ----
      #**********************************
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Group overlapping bacterial species ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, 
                                    "Start grouping bacteria"))
      
      blastSegments = 
        read_csv(paste0(sample, "/blastSegments.csv"), 
                 show_col_types = FALSE) %>% 
        select(name, LN, KC)
      blastOut = blastOut %>% left_join(blastSegments, 
                                        by = c("segmentId" = "name")) %>% 
        mutate(KC = KC * coverage * ident)
      
      # blastOut = blastOut %>% group_by(clusterId) %>% 
      #   filter(bit_score >= 0.90 * max(bit_score)) %>% ungroup()
      
      
      allBact = blastOut %>%
        # filter(geneId == "253") %>%
        select(segmentId, geneId, bit_score, coverage,
               accession, taxid, genus, species, extra, plasmid, KC, LN) %>%
        group_by(segmentId, geneId, accession) %>%
        filter(bit_score == max(bit_score)) %>% 
        dplyr::slice(1) %>% ungroup() %>%
        mutate(start = str_detect(segmentId, '_start$')) %>% 
        #If a blast result has identical scored for the top bact,
        #then add this score to each accession of that species because
        #the results are likely cropped by too many idential values
        group_by(segmentId) %>% 
        mutate(
          x = n_distinct(bit_score),
          accession = if_else(x == 1, as.character(taxid), accession)
        ) %>% 
        ungroup() 
      
      
      allBact = bind_rows(
        allBact %>% filter(x > 1),
        allBact %>% filter(x == 1) %>% 
          select(-extra, -accession) %>% distinct() %>% 
          left_join(allBact  %>% filter(x > 1) %>% 
                      select(geneId, accession, taxid) %>% distinct(), 
                    by = c("geneId", "taxid")) %>% 
          filter(!is.na(accession))
      ) %>% select(-x) %>% 
        #Add the path data
        left_join(pathData %>% 
                    select(pathId, segmentId, order, startOrientation, dist) %>% 
                    mutate(startOrientation = ifelse(str_detect(segmentId, "_start$"),
                                                     -1, startOrientation)) %>% 
                    distinct(), by = "segmentId") %>% 
        #Remove those not in a path
        filter(!is.na(order)) %>% rowwise() %>% 
        # mutate(pathScore = bit_score * coverage / (max(dist,0) * 0.001 + 1)) %>%
        mutate(pathScore = bit_score * coverage) %>%
        group_by(pathId, geneId, accession, startOrientation) %>% 
        filter(all(2:max(order) %in% order) | all(order == 1)) %>%
        group_by(geneId, accession, taxid, order, startOrientation) %>% 
        #Get the best paths per accession
        filter(pathScore == max(pathScore)) %>% dplyr::slice(1) %>% 
        group_by(geneId, accession, taxid, genus, species, plasmid) %>% 
        filter(any(order == 1)) %>%
        ungroup()
      
      #Make sure that taxid always has the same genus / species name and vice versa
      allBact = allBact %>% 
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
      
      allBact = map_df(genesDetected %>% filter(keep, !is.na(bactGroup)) %>% 
                         pull(geneId), function(geneId){
        
        myGene = allBact %>% filter(geneId %in% {{geneId}})
        
        bactSegments = myGene %>% group_by(segmentId, taxid) %>% 
          summarise(.groups = "drop") %>% 
          group_by(taxid) %>% summarise(x = list(segmentId)) 
        
        gfaMatrix = sapply(bactSegments$x, function(x){
          len = length(x)
          sapply(bactSegments$x, function(y){
            sum(x %in% y) / len
          })
        })
        
        rownames(gfaMatrix) = bactSegments$taxid
        colnames(gfaMatrix) = bactSegments$taxid
        
        myGene = myGene %>% 
          group_by(geneId, accession, taxid, genus, species, plasmid) %>% 
          summarise(n = n(), 
                    extension = sum(pathScore[!start]),
                    pathScore = sum(pathScore), 
                    KC = sum(KC), LN = sum(LN),
                    .groups = "drop") %>%
          group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
          arrange(desc(pathScore))

        myTaxid = myGene$taxid %>% unique()
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
        
        myGene %>% filter(taxid %in% myTaxid) %>% 
          group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
          slice(1)
        
      }) %>% mutate(depth = KC / LN) %>% 
        left_join(
          genesDetected %>% mutate(geneId = as.character(geneId)) %>% 
            select(geneId, gene, subtype, cover1), 
          by = "geneId"
        )
      
      
      #If extension is 0, it's not likely connected to a bact with genes
      # that have great extensions, even is the bact is present in both
      
      test = allBact %>% group_by(geneId) %>% 
        mutate(val = pathScore / max(pathScore),
               depth = KC / LN)  %>% 
        left_join(
          genesDetected %>% mutate(geneId = as.character(geneId)) %>% 
            select(geneId, gene, subtype, cover1)
        )
      

      
      geneMatrix = allBact %>%
        group_by(geneId, bact = taxid) %>%
        summarise(pathScore = round(max(pathScore),2), 
                  extension = round(max(extension),2), 
                  .groups = "drop") %>%
        left_join(genesDetected %>%
                    select(geneId, ARGgroup, keep) %>%
                    mutate(geneId = as.character(geneId)), by = "geneId") %>%
        filter(keep)
      
      # geneMatrix = allBact %>%
      #   mutate(bact = paste(genus, species)) %>%
      #   group_by(geneId, bact) %>%
      #   summarise(pathScore = round(max(pathScore),2), .groups = "drop") %>%
      #   left_join(genesDetected %>%
      #               select(geneId, ARGgroup, keep) %>%
      #               mutate(geneId = as.character(geneId)), by = "geneId") %>%
      #   filter(keep)
      
      geneMatrix = geneMatrix %>%
        select(bact, geneId, pathScore) %>%
        pivot_wider(bact,
                    names_from = "geneId", values_from = "pathScore",
                    values_fill = 0) %>%
        column_to_rownames("bact") %>% as.matrix()
      
      #Normalise the matrix per column
      normalise = function(x){
        (x - min(x)) / (max(x) - min(x))
      }
      geneMatrix = apply(geneMatrix, 2, normalise)
      
      # --- test 1
      myAdjustment = data.frame(
        geneId = colnames(geneMatrix)
      ) %>% left_join(
        genesDetected %>% select(geneId, cover1) %>%
          mutate(geneId = as.character(geneId)), by = "geneId")

      geneMatrix = t(t(geneMatrix) * myAdjustment$cover1)
      # --- end test 1
      
      #Calulcate the adjustment for each cell based on other scores in the row
      # the higher the ther scores (e.g. other gened detected), the more scaled up
      myAdjustment = (matrix(rowSums(geneMatrix), nrow = nrow(geneMatrix), 
                             ncol = ncol(geneMatrix)) - geneMatrix) 
      myAdjustment = myAdjustment / 
        matrix(apply(geneMatrix, 2, function(x) min(x[x > 0])), 
               nrow = nrow(geneMatrix), ncol = ncol(geneMatrix), byrow = T)
      myAdjustment[myAdjustment == 0] = 1
      
      # geneMatrix = geneMatrix + (myAdjustment / max(myAdjustment)) * 0.5
      # geneMatrix = (geneMatrix)^2 * myAdjustment
      geneMatrix = geneMatrix * myAdjustment
      
      # --- test 2
      # myScale = rowSums(geneMatrix)
      # myScale = myScale / min(myScale)
      # geneMatrix = geneMatrix * myScale
      
      # geneMatrix = geneMatrix - ((geneMatrix * 100) %% 1) / 100
      # --- end test 2
      
      #Filter to call for each gene one bacterial cluster
      myClusters = apply(geneMatrix, 2, function(x) {x == max(x)})
      myClusters = apply(myClusters, 1, function(x) {colnames(myClusters)[x]})
      myClusters = myClusters[sapply(myClusters, length) > 0]
      names(myClusters) = 1:length(myClusters)
      
      #Get the list of bact that belong to a cluster and their probability
      bactList = map_df(myClusters, function(group){
        group = geneMatrix[,group] %>% as.data.frame()
        group$prob = rowSums(group)
        group %>% filter(prob > 0) %>% mutate(prob = softmax(prob)) %>% 
          rownames_to_column("bact") %>% select(bact, prob) 
      }, .id = "bactGroup") %>% 
        mutate(bact = as.integer(bact), bactGroup = as.integer(bactGroup)) %>% 
        left_join(allBact %>% select(bact = taxid, genus, species) %>% 
                    distinct(), by = "bact") %>% 
        group_by(bactGroup) %>% 
        filter(prob > min(
          sort(unique(prob), decreasing = T)[1:(min(11, n_distinct(prob)))])
        ) %>% ungroup() %>% 
        arrange(bactGroup, desc(prob)) %>% ungroup() %>% 
        mutate(taxId = bact, 
               species = str_replace(species, "sp\\d+", "sp.")) %>% 
        left_join(bactGenomeSize %>% select(genus, species, size) %>% distinct(),
                  by = c("genus", 'species')) %>% 
        mutate(val = NA, runId = {{runId}},
               pipelineId = toProcess$pipelineId[sampleIndex])
      
      #Convert the clusters list to a dataframe format
      myClusters = map_df(myClusters, function(bact){
        data.frame(geneId = bact)
      }, .id = "bactGroup")
      
      #In some cases there is no data from blast about a gene, 
      # the bacgroup will become NA (warning suppressed)
      genesDetected = suppressWarnings(genesDetected %>% 
        select(-matches("bactGroup")) %>% 
        left_join(myClusters %>% 
                    mutate(geneId = as.integer(geneId)), by = "geneId") %>% 
        group_by(ARGgroup) %>% 
        mutate(bactGroup = as.integer(min(bactGroup, na.rm = T))) %>% 
        ungroup())
      
      #Calculate the RA
      bactList = bactList %>% 
        left_join(
              genesDetected %>% select(bactGroup, startDepth) %>%
                group_by(bactGroup) %>%
                summarise(estimatedAbundance = mean(startDepth), 
                          .groups = "drop"),
              by = "bactGroup") %>%
        mutate(estimatedAbundance = estimatedAbundance * size * 1e6 / inputfileBP,
               across(c(prob, estimatedAbundance), round, digits = 4),
               runId = {{runId}}) %>% select(-bact)
      
      
      # bact = result %>% group_by(from, bactGroup, genus, species) %>%
      #   summarise(val = max(val), .groups = "drop") %>%
      #   group_by(bactGroup) %>%
      #   mutate(prob = softmax(val, T, T)) %>%
      #   ungroup() %>%
      #   mutate(
      #     pipelineId = toProcess$pipelineId[sampleIndex]) %>%
      #   arrange(bactGroup, desc(prob)) %>% mutate(
      #     species = str_replace(species, "sp\\d+", "sp.")
      #   ) %>% left_join(bactGenomeSize %>% select(genus, species, level, size),
      #                   by = c("genus", "species")) %>%
      #   select(pipelineId, taxId = from, bactGroup, genus, species, prob, val, size, level) %>%
      #   left_join(
      #     genesDetected %>% select(bactGroup, startDepth) %>%
      #       group_by(bactGroup) %>%
      #       summarise(estimatedAbundance = mean(startDepth)),
      #     by = "bactGroup") %>%
      #   mutate(estimatedAbundance = estimatedAbundance * size * 1e6 / inputfileBP,
      #          across(c(prob, estimatedAbundance, val), round, digits = 4),
      #          runId = runId) %>% group_by(pipelineId, taxId) %>%
      #   filter(val == max(val)) %>% slice(1) %>% ungroup()
      
     
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, 
                                    "Finished grouping bacteria"))
      
      
      # ---- MAKE PREDICTIONS ----
      #***************************
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Predicting AMR ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, 
                                    "Start AMR prediction"))
      
      myInputs = genesDetected %>%
        mutate(id = paste(pipelineId,  bactGroup, sep = "_")) %>%
        group_by(pipelineId, bactGroup, gene) %>%
        mutate(val = startPerc * startDepth * cover1) %>% 
        filter(val == max(val)) %>% dplyr::slice(1) %>%
        ungroup() %>% 
        select(id, gene, cover1) %>%
        pivot_wider(gene, names_from = id, values_from = cover1)
      
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
      
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
      sqliteSetBusyHandler(myConn, 5000)
      
      #Add the ARGgroup to the detectedARG table
      q = dbExecute(myConn, 
                    "UPDATE detectedARG 
            SET runId = ?, ARGgroup = ?, bactGroup = ?
            WHERE pipelineId = ? AND geneId = ?", 
                    params = list(genesDetected$runId, genesDetected$ARGgroup, 
                                  genesDetected$bactGroup,
                                  genesDetected$pipelineId, genesDetected$geneId))
      x = q
      #Add the detected bacteria to the detectedBact table
      q = dbExecute(
        myConn, 
        sprintf("DELETE FROM detectedBact WHERE pipelineId = %i", bactList$pipelineId[1]))
      
      q = dbExecute(myConn, 
                    "INSERT INTO detectedBact (pipelineId, runId, taxId, bactGroup, 
            genus, species, prob, relAbundance, value) 
            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", 
                    params = bactList %>% 
                      select(pipelineId, runId, taxId, bactGroup, genus, species, 
                             prob, estimatedAbundance, val) %>% 
                      as.list() %>% unname())
      x = c(x, q)
      #Add the AMT predictions bacteria to the AMRprediction table
      q = dbExecute(
        myConn, 
        sprintf("DELETE FROM AMRprediction WHERE pipelineId = %i", bactList$pipelineId[1]))
      q = dbExecute(myConn, 
                    "INSERT INTO AMRprediction (pipelineId, runId, 
            bactGroup, antibioticId, resistance, value) 
            VALUES(?, ?, ?, ?, ?, ?)", 
                    params = predictions %>% 
                      select(pipelineId, runId, bactGroup, antibioticId,
                             resistance, value) %>% 
                      as.list() %>% unname())
      x = c(x, q)
      dbDisconnect(myConn)
      cat(myPipelineId, x, "\n")
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
      
      
    },
    finally = {
      
      #Write the logs to the database, irrespective of success
      newLogs$runId = runId
      newLogs$tool = "ARG_annotation.R"
      newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
      
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
      sqliteSetBusyHandler(myConn, 5000)
      
      q = dbExecute(
        myConn, 
        "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
        params = unname(as.list(newLogs)))
      
      #If pipeline successfully completed, note that in DB
      if(14 %in% newLogs$actionId){
        q = dbExecute(
          myConn, 
          "UPDATE pipeline SET statusCode = 5, inputFileBP = ?,
           statusMessage = 'Finished pipeline', modifiedTimestamp = ?
           WHERE pipelineId = ?",
          params = list(inputfileBP, as.character(Sys.time()),myPipelineId))
      }
      
      dbDisconnect(myConn)
      
    })
    
  }
  # )
  # rm(cl)
}

