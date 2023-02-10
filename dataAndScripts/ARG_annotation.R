#**************************
# ---- ARG ANNOTATION ----
#*************************

suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(RSQLite))
# suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(rmarkdown))
# suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

args = commandArgs(trailingOnly = TRUE) 
baseFolder = formatPath(args[[1]], endWithSlash = T)
database = args[[2]]
runId = as.integer(args[[3]])
verbose = abs(as.integer(args[[4]]))
pipelineIds = str_trim(unlist(strsplit(args[[5]], ",")))
generateReport = as.logical(args[[6]])
maxCPU = as.integer(args[[7]])

# forceRedo = as.logical(as.logical(args[[7]]))
forceRedo = T
minBlastLength = 250

#Load the ARG and the sample list to process
myConn = dbConnect(SQLite(), database,  synchronous = NULL)

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

bactGenomeSize = 
  read.csv(sprintf("%sdataAndScripts/bactGenomeSize.csv",baseFolder))

dbDisconnect(myConn)

#Check if the temp folders still exist
toProcess$exists = dir.exists(toProcess$outputFolder)
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

blast_readOutput = function(file, separate = T, includeIssues = F, verbose = 1){
  
  outfmt = "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart send bitscore score length pident nident qcovs qcovhsp"
  
  #Load the output csv (zipped)
  blastOut = read.table(file, sep = "\t", quote = "", comment.char = "")
  colnames(blastOut) = strsplit(outfmt, " ")[[1]][-1]
  
  #split multiple matches
  blastOut = blastOut %>% 
    mutate(x = str_count(staxids, ";"), y = str_count(sallacc, ";"), 
           z = x == y | x == 0, 
           plasmid = str_detect(salltitles, "lasmid |ntegron |ransposon |ransposase ")) 

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
  # registerDoParallel(cores=min(maxCPU, length(sampleIndex)))
  result = foreach(sampleIndex = sampleIndex, 
                   packages = c("gfaTools","tidyverse","RSQLite")) %do% {
                     
    newLogs = data.frame(
      timeStamp = as.integer(Sys.time()), 
      actionId = 1, actionName = "Start Annotation")
                     
    #Process each sample
    tryCatch({

      sample = toProcess$outputFolder[sampleIndex]
      myPipelineId = toProcess$pipelineId[sampleIndex]
      sampleName = toProcess$name
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          sprintf("Start Annotation for pipeline ID %i\n", 
                                  myPipelineId))}
      
      inputfileBP = grep("sequences length",
                         readLines(sprintf("%s/metacherchant_logs/log", sample)),
                         value = T, fixed = T) %>%
        str_extract("([\\d'])+(?=\\s\\()") %>%
        str_replace_all("'", "") %>% as.numeric()
      
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
        blast_readOutput) %>% bind_rows() 
      
      if(nrow(blastOut) == 0) {
        allBact = data.frame()
        stop()
      }
      
      expandBlast = list.files(sample, full.names = T, pattern = "expand_\\d+.csv.gz")
      
      if(length(expandBlast) > 0){
        
        #Get the new results
        blastOut2 = lapply(
          list.files(sample, full.names = T,
                     pattern = "expand_\\d+.csv.gz"),
          blast_readOutput) %>% 
          bind_rows()
        
        blastOut = bind_rows(
          blastOut %>% filter(!qseqid %in% unique(blastOut2$qseqid)),
          blastOut2
        )
        
        rm(blastOut2)
      }
      
      #Extract data we need + transform
      blastOut = blastOut %>% 
        select(query_title = qseqid, taxid = staxids, accession = sallacc, 
               bact = sscinames, plasmid, bit_score = bitscore, score, 
               identity = nident, query_len = qlen, query_from = qstart, query_to = qend, 
               hit_from = sstart, hit_to = send, align_len = length) %>% 
        mutate(bact = str_remove_all(bact, "[^\\w\\s]")) %>% 
        extract(bact, c("genus", "species", "extra"), 
                regex = "(\\w+)\\s+(\\w+)($|\\s+.*)") %>% 
        #Sp. will be pasted with taxid to make it unique
        mutate(species = ifelse(species == "sp", paste0(species, taxid), species)) %>% 
        filter(!genus %in% c("uncultured", "Uncultured", "mixed"), !is.na(species)) %>%
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
          coverage = min(align_len / query_len, 1), #In case of gap can be > 1
        ) %>% ungroup() 
      
      #Fix some known issues with naming
      blastOut[blastOut$genus == "Enterobacter" & 
                 blastOut$species == "aerogenes","genus"] = "Klebsiella"
      
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                    "Finished loading BLASTn output"))
      
      #Check
      blastSegments = 
        read.csv(paste0(sample, "/blastSegments.csv")) %>% 
        select(name, LN, KC) %>% 
        mutate(
          start = ifelse(str_detect(name, "_start$"), T, F),
          geneId = as.integer(str_extract(name, "\\d+"))) %>% 
        #Correct KC where depth is an outlier
        group_by(geneId) %>% mutate(
          depth = KC / LN,
          z = (depth - mean(depth)) / sd(depth),
          depth = ifelse(z > 1.96,  sum(KC[z <= 1.96 & LN >= 250]) / 
                           sum(LN[z <= 1.96 & LN >= 250]), depth),
          depth = ifelse(is.na(depth) | is.nan(depth), KC / LN, depth),
          KC = depth * LN
        ) %>% ungroup() %>% select(-z, -depth)
      
      blastOut = blastOut %>% 
        left_join(blastSegments, by = c("segmentId" = "name")) %>% 
        mutate(KC = KC * coverage * ident) #Lower coverage is lower KC score
      

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
      
      
      #Interpolate bit_score FOR non-blasted < 250 pieces
      #Calculate the score to bit_score conversion factor for each gene
      extraBits = blastOut %>% filter(ident == 1, coverage == 1) %>% 
        group_by(geneId) %>% 
        summarise(bitConst = mean(bit_score / score), .groups = "drop")
      
      extraBits = pathData %>% 
        mutate(LN = ifelse((LN < minBlastLength & dist >= 0 ) | (dist < 0 & LN < 75),
                           LN, 0)) %>% 
        group_by(geneId, orientation, pathId) %>% 
        transmute(order, score = cumsum(LN) - 30*(1:n()),
                  score = ifelse(score < 0 , 0, score)) %>% 
        ungroup() %>% 
        left_join(extraBits, by = "geneId") %>% 
        mutate(bitScore = bitConst * score)
      
      #Only work with segments >= 250 in paths
      pathData = pathData %>% 
        filter(LN >= minBlastLength | (dist < 0 & LN > 74)) %>% 
        group_by(geneId, pathId) %>% 
        arrange(geneId, pathId, desc(order)) %>% 
        mutate(oldOrder = order, order = n():1) %>% 
        ungroup()
      
      #Also add the fragment data to the paths
      fragments = gfa_read(paste0(sample, "/fragmentGFA.gfa"))
      
      if(nrow(fragments$segments) > 0){
        #Fragments were merged, so treated as a non-start piece
        fragments = fragments$segments %>% 
          select(segmentId = name, LN, KC, geneId) %>% 
          group_by(geneId) %>% 
          mutate(pathId = 0, orientation = 0:(n()-1),
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
      #Grab the detected ARG from the previous step
      myConn = dbConnect(SQLite(), database,synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      
      genesDetected = tbl(myConn, "detectedARG") %>% 
        filter(pipelineId == myPipelineId) %>% as.data.frame() %>% 
        mutate(cover = ifelse(type == 'noFragments', cover1, cover2))
      dbDisconnect(myConn)
      
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

      #Add the group (membership of subgraphs)
      blastOut = blastOut %>% left_join(
        pathData %>% select(segmentId, group) %>% distinct(), by = "segmentId") %>% 
        mutate(geneId = as.integer(geneId))

      #Find duplicate genes if they lie on same pos in genome
      toCheck = blastOut %>% 
        group_by(geneId, accession, group) %>% 
        filter(any(start)) %>% 
        mutate(sMin = ifelse(hit_from < hit_to, hit_from, hit_to),
               sMax = ifelse(hit_from >= hit_to, hit_from, hit_to)) %>% 
        ungroup()
      
      toCheck = bind_rows(
        toCheck %>% filter(start) %>% 
          select(accession, geneId, pos = sMin, clusterId) %>% 
          group_by(accession, pos) %>% filter(n() > 1) %>% 
          ungroup(),
        toCheck %>% filter(start) %>% 
          select(accession, geneId, pos = sMax, clusterId) %>% 
          group_by(accession, pos) %>% filter(n() > 1) %>% 
          ungroup()
      ) %>%
        group_by(accession) %>% 
        mutate(n = n()) %>% ungroup() %>% 
        select(accession, geneId, pos, n, clusterId) %>% 
        left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>% 
        arrange(desc(n), accession, pos)
      
      duplicates = data.frame(geneId = integer(), cluster = double())
      
      while(nrow(toCheck) > 0){
        nextCheck = toCheck %>% 
          filter(accession == accession[1]) 
        
        lastgroup = nextCheck$pos
        nextCheck = nextCheck %>% 
          group_by(geneId) %>% mutate(group = min(pos)) %>% 
          group_by(pos) %>% mutate(group = min(group)) %>% ungroup()
        
        while(any(lastgroup != nextCheck$group)){
          lastgroup = nextCheck$group
          nextCheck = nextCheck %>% 
            group_by(geneId) %>% mutate(group = min(group)) %>% 
            group_by(pos) %>% mutate(group = min(group)) %>% ungroup()
        }
        
        nextCheck = nextCheck %>% select(geneId, group) %>% distinct() %>% 
          left_join(duplicates, by = "geneId") %>% 
          group_by(group) %>% 
          mutate(cluster = ifelse(any(!is.na(cluster)), min(cluster, na.rm = T), group[1])) %>% 
          ungroup() %>% select(-group)
        
        duplicates = bind_rows(duplicates, nextCheck) %>% distinct()
        
        toCheck = toCheck %>% group_by(accession) %>% 
          filter(!all(geneId %in% duplicates$geneId)) %>% ungroup() %>% 
          arrange(desc(n), accession, pos)
      }
      
      duplicates = duplicates %>% 
        transmute(geneId, duplicated = as.factor(cluster) %>% as.integer()) 
      
      #Add the clusters to detected ARG and fill in blanks
      genesDetected = genesDetected %>% 
        left_join(duplicates, by = "geneId") 
      
      if(any(is.na(genesDetected$duplicated))){
        newIdx = max(genesDetected$duplicated, 0, na.rm = T) + 1
        genesDetected$duplicated[is.na(genesDetected$duplicated)] = 
          newIdx:(sum(is.na(genesDetected$duplicated))+newIdx-1)
      }
      
      blastOut$start[blastOut$start & ! blastOut$segmentId %in% 
                       segmentsOfInterest$name[segmentsOfInterest$type != "fragmentsOnly"]] = F
      
      
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
        filter(all(order != 1))
      if(nrow(addSeedHit) > 0){
        addSeedHit = addSeedHit %>%
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
          left_join(addSeedHit %>% 
                      select(geneId, accession, hit_from, hit_to, group),
                    by = "geneId") %>% 
          left_join(
            blastOut %>% 
              select(accession, taxid, genus, species, extra, 
                     subspecies, plasmid) %>% 
              distinct(), by  = "accession")
        
        #Add the new start 'matches' to the blast output
        blastOut = bind_rows(blastOut, addSeedHit)
      }
      #GENOME DIST START CHECK
      
      #Order the min - max positino of match in target genome
      tempVar = blastOut %>% 
        group_by(geneId, accession, group) %>% 
        filter(any(start)) %>% 
        mutate(sMin = ifelse(hit_from < hit_to, hit_from, hit_to),
               sMax = ifelse(hit_from >= hit_to, hit_from, hit_to))
      
      #Get the locations of all start matches in same genome
      tempVar2 = tempVar %>% filter(start) %>% 
        select(geneId, accession, group, sMin) %>% 
        mutate(startGroup = 1:n()) %>% 
        pivot_wider(names_from = startGroup, values_from = sMin) %>% 
        ungroup()
      
      #Add them to the data frame
      tempVar  = tempVar %>% ungroup() %>% 
        left_join(tempVar2, by = c("geneId", "accession", "group"))
      
      #For each segment, pick the closest start segment for comparison
      tempVar$startGroup = 
        apply(tempVar %>% select(sMin, matches(c("1", "2", "3"))), 1, function(x){
          x = abs(x[-1] - x[1])
          which(x == min(x, na.rm = T))
        })
      
      #Only keep a segment if it is the correct distance from the start segment
      tempVar = tempVar %>% 
        group_by(geneId, accession, group, startGroup) %>% 
        filter(any(start)) %>% 
        mutate(
          after = ifelse(sMin > sMin[start], T, F),
          dist = ifelse(after, sMin - sMax[start], sMin[start] - sMax) + 29) %>%
        ungroup() %>% 
        left_join(pathData %>% select(segmentId, dist1 = dist) %>% 
                    distinct(), by = "segmentId") %>% 
        rowwise() %>% 
        #Repeats can make segments shorter because loop is cut off, so provide some slack
        mutate(correct = between(dist, dist1 -500, dist1 + 500) | start) %>% 
        filter(correct) %>% 
        group_by(geneId, accession, group, startGroup) %>% 
        mutate(bestMatch = sum(bit_score)) %>% 
        group_by(geneId, accession, group) %>% 
        filter(bestMatch == max(bestMatch)) %>% 
        group_by(segmentId, accession) %>% 
        filter(bit_score == max(bit_score)) %>% 
        slice(1) %>% ungroup()
      
      wrongLocation = blastOut %>% filter(!segmentId %in% tempVar$segmentId)
      
      blastOut = bind_rows(
        tempVar %>% select(matches(colnames(blastOut))),
        
        blastOut %>% 
          filter(geneId %in% (genesDetected %>% 
                                filter(!geneId %in% segmentsOfInterest$geneId) %>% 
                                pull(geneId))) %>% 
          group_by(segmentId, accession) %>% filter(bit_score == max(bit_score,0)) %>% 
          slice(1) %>% ungroup(),
        
        blastOut %>% 
          filter(segmentId %in% 
                   fragments$segmentId[str_detect(fragments$segmentId, "_start$")]) %>% 
          group_by(segmentId, accession) %>% filter(bit_score == max(bit_score,0)) %>% 
          slice(1) %>% ungroup()
      )
      
      blastOut$start[blastOut$segmentId %in% 
                       fragments$segmentId[str_detect(fragments$segmentId, "_start$")]] = T
      
      rm(tempVar, tempVar2)
      
      allBact = blastOut %>%
        select(segmentId, geneId, bit_score, coverage,
               accession, taxid, genus, species, extra, plasmid, KC, LN, start) %>%
        
        #Only keep best match per segmentId per accession
        group_by(segmentId, geneId, accession) %>%
        filter(bit_score == max(bit_score)) %>% 
        dplyr::slice(1) %>% ungroup() %>%
        
        #Adjust path score based on coverage (???)
        mutate(pathScore = bit_score * coverage) %>%
        
        #Add the path data and only keep consecutive segments (ignore < 250)
        left_join(pathData %>% 
                    select(pathId, segmentId, order, oldOrder, orientation) %>% 
                    mutate(orientation = ifelse(str_detect(segmentId, "_start$"),
                                                -1, orientation)) %>% 
                    distinct(), by = "segmentId") %>% 
        filter(!is.na(order))  %>% 
        arrange(geneId, pathId, accession, plasmid, order) %>% 
        group_by(geneId, pathId, accession, plasmid) %>% 
        filter(order == 1:n() | orientation == -1) %>% 
        
        #Calculate the total path score
        group_by(geneId, accession, plasmid, orientation, pathId) %>%
        mutate(fullPath = sum(pathScore)) %>% ungroup()
      
      #Add the extra bits to each fullPath score
      tempVar = allBact %>% filter(orientation != -1) %>% 
        select(geneId, pathId, oldOrder, accession) %>% distinct() %>% 
        left_join(extraBits, by = c("geneId", "pathId", "oldOrder" = "order")) %>% 
        select(geneId, accession, pathId, oldOrder, extraBits = bitScore)
      
      allBact = allBact %>% left_join(
        tempVar, by = c("geneId", "pathId", "oldOrder", "accession")
      ) %>% group_by(geneId, accession, pathId, orientation) %>% 
        mutate(extraBits = max(extraBits, 0, na.rm = T),
               fullPath = fullPath + extraBits,
               extraBits = c(rep(0, n() - 1), extraBits[n()])) %>% ungroup()
      
      allBact = allBact %>%
        
        #Pick best path for each orientation
        group_by(geneId, accession, plasmid, orientation) %>%
        filter(pathId == pathId[fullPath == max(fullPath)][1]) %>%
        
        #Sum the 3 orientations to get full path score
        group_by(geneId, accession, plasmid) %>%
        mutate(fullPath = sum(pathScore) + sum(extraBits)) %>%
        group_by(geneId, accession) %>%
        filter(plasmid == plasmid[fullPath == max(fullPath)][1]) %>%
        ungroup()
      
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
        ungroup() 
      
      
      allBact = allBact %>% 
        group_by(geneId, accession, taxid, genus, species, plasmid, LN, KC) %>% 
        summarise(extension = sum(pathScore[!start]) + sum(extraBits),
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
      
      
      # Bact on the same path represent the same species
      allBact = allBact %>% group_by(geneId, path0, path1) %>% 
        mutate(pathGroup = cur_group_id())
      allBact = allBact %>% group_by(geneId, pathGroup) %>% 
        mutate(maxGroupScore = max(fullPath)) %>% 
        ungroup()
      
      # Add plasmid prob per gene
      allBact = allBact %>% 
        group_by(geneId) %>% mutate(top = fullPath  / maxGroupScore) %>% ungroup() %>% 
        left_join(
          allBact %>% group_by(geneId, plasmid) %>% 
            summarise(val = max(fullPath), .groups = "drop") %>% 
            group_by(geneId) %>% summarise(
              probPlas = max(0, val[plasmid]) / sum(val), .groups = "drop"
            ), by  = "geneId")
      
      allBact = allBact %>% 
        mutate(pipelineId = myPipelineId, runId = runId) %>% 
        select(-c(gene:type, maxGroupScore:probPlas))
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, 
                                    "Finished grouping"))
      
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
            SET runId = ?, ARGgroup = ?
            WHERE pipelineId = ? AND geneId = ?", 
                    params = list(genesDetected$runId, genesDetected$duplicated,
                                  genesDetected$pipelineId, genesDetected$geneId))
      q = dbClearResult(q)
      
      #Add new bacterial taxa and strains
      newStrains = allBact %>% select(accession, taxid, plasmid) %>% 
        filter(!accession %in% (tbl(myConn, "bactStrains") %>% pull(accession))) %>% 
        mutate(taxid = as.integer(taxid)) %>% distinct()
      
      if(nrow(newStrains) > 0){
        q = dbSendStatement(myConn, 
            "INSERT INTO bactStrains (accession, taxid, plasmid) 
            VALUES(?, ?, ?)", 
            params = newStrains %>% as.list() %>% unname())
        q = dbClearResult(q)
      }
      
      newTaxa = allBact %>% select(taxid, genus, species) %>% distinct() %>% 
        mutate(taxid = as.integer(taxid)) %>% 
        filter(!taxid %in% (tbl(myConn, "bactTaxa") %>% pull(taxid))) 
      
      if(nrow(newTaxa) > 0){
        q = dbSendStatement(myConn, 
            "INSERT INTO bactTaxa (taxid, genus, species) 
            VALUES(?, ?, ?)", 
            params = newTaxa %>% as.list() %>% unname())
        q = dbClearResult(q)
      }

      #Add the new annotation data
      q = dbSendStatement(
        myConn, 
        sprintf("DELETE FROM annotation WHERE pipelineId = %i", 
                myPipelineId))
      q = dbClearResult(q)
      
      q = dbSendStatement(myConn, 
            "INSERT INTO annotation (pipelineId, runId, geneId, accession,
            LN, KC, extension, fullPath, extraBits, path0, maxOrder0, path1, 
            maxOrder1, pathGroup) 
            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", 
            params = allBact %>% 
              mutate(across(c("KC", "extension", "fullPath", "extraBits"), round, digits = 2)) %>% 
              select(pipelineId, runId, geneId, accession, LN, KC, extension, 
                     fullPath, extraBits, path0, maxOrder0, path1, maxOrder1, 
                     pathGroup) %>% as.list() %>% unname())
      q = dbClearResult(q)
      
      #Write output file to be used as alternative input for Explorer App
      print(generateReport)
      if(generateReport){
        td = tempdir()
        dir.create(sprintf("%s/results_id%i", td, myPipelineId), F)
        
        file1 = tbl(myConn, "pipeline") %>%
          filter(pipelineId == myPipelineId) %>% 
          collect()
        write_csv(file1, sprintf("%s/results_id%i/info.csv", td, myPipelineId))
        
        
        file2 = tbl(myConn, "detectedARG") %>%
          filter(pipelineId == myPipelineId) %>% 
          # select(-gene, -subtype) %>% 
          collect() %>% 
          left_join(ARG %>% select(-type), by = "geneId") %>% 
          mutate(
            cover = ifelse(type == "noFragments", cover1, cover2),
            across(c(cover, startPerc), function(x) round(x * 100, 2)),
            startDepth = round(startDepth, 2)) %>% 
          arrange(desc(cover), desc(startDepth))
        write_csv(file2, sprintf("%s/results_id%i/detectedARG.csv", td,myPipelineId))
        
        file3 = tbl(myConn, "annotation") %>%
          filter(pipelineId == myPipelineId) %>%
          left_join(tbl(myConn, "bactStrains"), by = "accession") %>%
          left_join(tbl(myConn, "bactTaxa"), by = "taxid") %>% 
          collect() %>% 
          group_by(geneId) %>% mutate(
            top = fullPath / max(fullPath),
            depth = KC / LN
          ) %>% ungroup() %>% left_join(
            file2 %>% select(geneId, gene, subtype, cover, type), 
            by = "geneId") %>% 
          mutate(plasmid = as.logical(plasmid))
        write_csv(file3, sprintf("%s/results_id%i/annotation.csv", td, myPipelineId))
        
        system(sprintf("tar -czvf %s --directory %s %s",
               sprintf("%s/results_id%i.tar.gz", sample, myPipelineId),td,
               sprintf("results_id%i", myPipelineId)), intern = F)
      }
      
      dbDisconnect(myConn)
      
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
      
    })
    
  }
}

