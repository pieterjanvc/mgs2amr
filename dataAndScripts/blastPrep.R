#!/usr/bin/env Rscript

#**********************
# ---- Blast prep ----
#*********************
require(igraph)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))


# ---- Inputs ----
#*****************
args = commandArgs(trailingOnly = TRUE) #arguments specified in the shell file

baseFolder = formatPath(args[[1]], endWithSlash = T)
database = args[[2]]
tempFolder = formatPath(args[[3]], endWithSlash = T)
tempName = args[[4]]
verbose = abs(as.integer(args[[5]]))
runId = as.integer(args[[6]])
pipelineId = as.integer(args[[7]])
keepAllMetacherchantData = as.logical(args[[8]])
maxPathDist = as.integer(args[[9]]) #Distance from ARG to crop the GFA file (reduces blast search)
minBlastLength = as.integer(args[[10]]) #Min segment length to submit to blast
trimLength = as.integer(args[[11]]) #Loose segments smaller than this will be cut from thr GFA
clusterIdentidy  = as.numeric(args[[12]]) #The cluster identity percent used in usearch
forceRedo = as.logical(args[[13]]) #If parts of the code have successfully run before a crash, do not repeat unless forceRedo = T
maxStep = as.integer(args[[14]]) #Which parts of the script to run? If NA all is run

maxCPU = 4
maxStep = ifelse(maxStep == 0, 5, maxStep)
maxStartGap = 800
minCover = 0.25 #The minimum coverage to conksider a ARG being present
maxPathIter = 8000

# database = sprintf("%sdataAndScripts/mgs2amr.db", baseFolder)
# database = "/mnt/meta2amrData/pipelineTest/after1200/mgs2amr.db"

#Generate a list out of the settings file
settings = readLines(paste0(baseFolder, "settings.txt"))
settings = settings[str_detect(settings,"^\\s*[^=#]+=.*$")]
settings = str_match(settings, "\\s*([^=\\s]+)\\s*=\\s*(.*)")
settings = setNames(str_trim(settings[,3]), settings[,2])

#Check if pigz is available instead of gzip for faster zipping
zipMethod = ifelse(length(suppressWarnings(
  system("command -v pigz", intern = T))) == 0,
  "gzip", "pigz")

tempFolder = formatPath(paste0(tempFolder, tempName), endWithSlash = T)
logPath = sprintf("%s%s_log.csv", tempFolder,tempName)

#Check the log file to see if there was a previous run of the code
myConn = dbConnect(SQLite(), database, synchronous = NULL)
sqliteSetBusyHandler(myConn, 30000)
logs = dbGetQuery(
  myConn,
  paste("SELECT * FROM logs WHERE runId IN",
        "(SELECT runId FROM scriptUse WHERE pipelineId == ?)",
        "AND tool = 'blastPrep.R'"),
  params = pipelineId)
dbDisconnect(myConn)
newLogs = data.frame(timeStamp = as.integer(Sys.time()),
                     actionId = 1, actionName = "Start BLAST prep")

switchLinks = function(links){
  
  if(nrow(links) == 0) return(links)
  
  links %>% mutate(
    fromOrient = ifelse(.$toOrient == "+", "-", "+"),
    toOrient = ifelse(.$fromOrient == "+", "-", "+"),
    from = .$to,
    to = .$from
  )
  
}

#! Added only for benchmarking
myStats = list(info = data.frame(
  pipelineId = pipelineId, runId = runId, 
  timeStamp = as.integer(Sys.time()), tempFolder = tempFolder))
#!

tryCatch({
  
  #Check if data exists
  if(length(list.files(tempFolder, pattern = "masterGFA.db")) == 0){
    cat("\nMGS2AMR cannot locate the data needed associated with this pipelineId",
        "\n Expected temp folder:", tempFolder, "\n\nPipeline halted\n")
    
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 23,
                                  "Data not found, halt pipeline"))
    maxStep = 0
    genesDetected = data.frame()
  }

  if(maxStep > 0){

    # ---- 1. Clean up files and folders ----
    #****************************************
    if(nrow(logs %>% filter(actionId %in% c(2, 6))) > 0){

      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "Skip MetaCherchant cleanup, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2,
                                    "Skip MetaCherchant cleanup, already done"))

      #Takes long time to load, only do if next step is not completed (not needed afterwards)
      # if(nrow(logs %>% filter(actionId %in% c(9,11))) == 0 | forceRedo){
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                            "Load master GFA file for processing ... ")}
        
      #! Get the database size 
      myStats$DBzipped = file.info(paste0(tempFolder, "masterGFA.db.gz"))$size
      #!
    
      system(sprintf("%s -d %s", zipMethod, paste0(tempFolder, "masterGFA.db.gz")))
      myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
      gfa = list(
        segments = dbReadTable(myConn, "segments"),
        links = dbReadTable(myConn, "links")
      )
      dbDisconnect(myConn)
      
      #! Get the masterGFA info 
      myStats$masterGFA = list(
        segments = gfa$segments %>% 
          mutate(start = str_detect(name, "_start$")) %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start), avgSLN = mean(LN[start]), avgSKC = mean(KC[start]),
            maxSLN  = max(LN[start]), .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId),
        links = gfa$links %>% group_by(geneId) %>% summarise(
          nLinks = n(),.groups = "drop"
        ) %>% mutate(pipelineId = pipelineId)
      )
      #!
      
      if(verbose > 0){cat("done\n")}
      

    } else { #Not processed yet ...

      # ---- Merge all GFA files----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Filter and merge MetaCherchant output ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, "Start filtering and merging MetaCherchant output"))

      myFiles = list.files(sprintf("%s", tempFolder),
                           ".gfa", recursive = T, full.names = T)
      myFiles = myFiles[!str_detect(myFiles, "masterGFA")]

      #This process can be done in parallel so speed things up
      #Read all GFA files
      registerDoParallel(cores=maxCPU)
      gfa = foreach(myFile = myFiles) %dopar% {

        geneId = str_extract(myFile, "\\d+(?=/graph.gfa)")
        myGFA = gfa_fixMetacherchant(myFile)

        #Check if the file is not empty
        if(nrow(myGFA$segments) > 0){

          myGFA$segments$geneId = geneId

          if(nrow(myGFA$links) > 0){
            myGFA$links$geneId = geneId
          }

          return(myGFA)

        } else {
          return(NULL)
        }

      }

      #Merge all the GFAs that we need for further analysis
      myFile = sapply(gfa, is.null)
      notUsed = gfa[myFile] %>% unlist()
      gfa = gfa[!myFile]
      gfa = list(
        segments = bind_rows(sapply(gfa, "[", 1)),
        links = bind_rows(sapply(gfa, "[", 2))
      )
      rm(myFile)

      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4,
                                    "Finished filtering and merging MetaCherchant output"))
      
      #! Get the masterGFA info 
      myStats$masterGFA = list(
        segments = gfa$segments %>% 
          mutate(start = str_detect(name, "_start$")) %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start), avgSLN = mean(LN[start]), avgSKC = mean(KC[start]),
            .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId),
        links = gfa$links %>% group_by(geneId) %>% summarise(
          nLinks = n(),.groups = "drop"
        ) %>% mutate(pipelineId = pipelineId)
      )
      #!

      # ---- Update the master GFA file and zip it to save space ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Write master GFA ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Start writing master GFA"))

      myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
      dbWriteTable(myConn, "segments", gfa$segments, overwrite = T)
      dbWriteTable(myConn, "links", gfa$links, overwrite = T)
      dbDisconnect(myConn)

      #Remove Metacherchant Data if set
      if(keepAllMetacherchantData == F & (nrow(logs %>% filter(actionId == 4)) == 0)){
        system(sprintf("ls -d %s/* | grep -P \"/\\d+$\" | xargs rm -R", tempFolder))
      }

      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6,
                                    "Finished writing master GFA to zip"))

    }
  }


  if(maxStep > 1){

    # ---- 2. Detect important ARG ----
    #**********************************
    if(nrow(logs %>% filter(actionId %in% c(9, 11))) > 0 & !forceRedo){

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Skip ARG detection, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Skip ARG detection, already done"))

      genesDetected = read.csv(paste0(tempFolder, "genesDetected/genesDetected.csv"))
      fragments = gfa_read(paste0(tempFolder, "fragmentGFA.gfa"))
      segmentsOfInterest = read_csv(paste0(tempFolder, "segmentsOfInterest.csv"))


    } else {

      # ---- Clean up and merge start segments ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Recover ARG seed sequences ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Start recovering ARG seed sequences"))

      #Split the file in groups to merge start segments (multiple can be done at once)
      myGroups = gfa$links %>% group_by(geneId) %>% summarise(n = n())

      myGroup = 1
      curSum = 0
      myResult = rep(0, nrow(myGroups))
      for(i in 1:nrow(myGroups)){

        curSum = curSum + myGroups$n[i]

        if(curSum >= 200000){
          curSum = 0
          myGroup = myGroup + 1
        }

        myResult[i] = myGroup
      }

      myGroups$group = myResult

      registerDoParallel(cores=min(maxCPU, myGroup, na.rm = T))
      
      #Read all GFA files
      gfa = foreach(myGroup = unique(myGroups$group)) %dopar% {

        myGenes = myGroups$geneId[myGroups$group == myGroup]

        myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
        gfaGroup = list(
          segments = dbGetQuery(
            myConn,
            sprintf('SELECT * FROM segments WHERE geneId IN (%s)',
                    paste(myGenes, collapse = ","))),
          links = dbGetQuery(
            myConn,
            sprintf('SELECT * FROM links WHERE geneId IN (%s)',
                    paste(myGenes, collapse = ","))))
        dbDisconnect(myConn)

        #Add the geneId to the segment/link names
        gfaGroup$segments = gfaGroup$segments %>%
          mutate(name = paste0(geneId, "_", name))
        gfaGroup$links = gfaGroup$links %>%
          mutate(
            from = paste0(geneId, "_", from),
            to = paste0(geneId, "_", to)
          )

        #Cut out small appendages from graphs to make joining easier
        gfaGroup = gfa_trimLooseEnds(
          gfaGroup, 75, keepRemoving = F,
          exclude = gfaGroup$segments$name[
            str_detect(gfaGroup$segments$name, "_start$")])

        #Remove single splits
        gfaGroup = gfa_filterSegments(
          gfaGroup,
          gfa_getOneSplit(gfaGroup,
                          except = gfaGroup$segments$name[
                            str_detect(gfaGroup$segments$name, "_start$")]),
          action = "remove"
        )

        #Merge the start segments
        mergeStartSegments(
          list(
            segments = gfaGroup$segments %>% filter(geneId %in% myGenes) ,
            links = gfaGroup$links %>% filter(geneId %in% myGenes)
          ), maxGap = maxStartGap)

      }

      #Merge the results
      gfa = list(
        segments = bind_rows(sapply(gfa, "[", 1)),
        links = bind_rows(sapply(gfa, "[", 2))
      )

      #Add the geneId back to the links and unitig names
      gfa$links$geneId = str_extract(gfa$links$from, "^\\d+")
      # gfa$segments$start[is.na(gfa$segments$start)] = ifelse(
      #   str_detect(gfa$segments$name[is.na(gfa$segments$start)], "_start$"),
      #   1,0
      # )

      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8,
                                    "Finished recovering ARG seed sequences"))
      
      #! Get the updates after start seg merge 
      myStats$mergedStart = gfa$segments %>% 
        filter(start > 0) %>% 
        group_by(geneId) %>% summarise(
          nStart = sum(start > 0), 
          avgSLN = mean(LN), avgSKC = mean(KC), maxSLN  = max(LN),
          .groups = "drop"
        ) %>% mutate(pipelineId = pipelineId)
      #!


      # ---- Estimate ARG presence by coverage and fragmentation----
      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Detect ARG in the data ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, "Start detecting ARG"))

      # ---- Detect the type of graph (fragmented) ----
      #************************************************
      
      #Get the membership of each component (subgraph id)
      membership = igraph::graph_from_data_frame(data.frame(
        from = gfa$links$from,
        to = gfa$links$to
      ), directed = F)
      
      membership = igraph::components(membership)$membership
      membership = data.frame(
        name = names(membership),
        group = membership,
        row.names = NULL
      )
      
      membership = bind_rows(
        membership,
        gfa$segments %>% select(name) %>% 
          filter(!name %in% c(membership$name)) %>% 
          mutate(group = (max(membership$group)+1):(max(membership$group)+n()))
      )
      
      gfa$segments = gfa$segments %>% select(-matches("group"))
      gfa$segments = gfa$segments %>% left_join(membership, by = "name") 
      
      #Get all start segments in 'from' column
      startConn = gfa$links %>% mutate(
        s1 = str_detect(from, "_start$"), s2 = str_detect(to, "_start$")
      )
      
      startConn[startConn$s2 & !startConn$s1,] = switchLinks(
        startConn[startConn$s2 & !startConn$s1,]
      )
      
      startConn = bind_rows(
        startConn %>% filter(s1 | s2) ,
        startConn %>% filter(s1, s2) %>% switchLinks()
      ) 
      
      #Check if the connections have a next connection
      x = gfa$links %>% 
        filter(from %in% startConn$to | to %in% startConn$to)
      x[!x$from %in% startConn$to,] = switchLinks(x[!x$from %in% startConn$to,])
      x = bind_rows(
        x, x %>% filter(str_detect(to, "_start$")) %>% switchLinks()
      ) 
      
      startConn = startConn %>% left_join(
        x %>% select(from, fromOrient) %>% distinct() %>% mutate(second = T), 
        by = c("to" = "from", "toOrient" = "fromOrient")) %>% 
        mutate(second = replace_na(second, F))
      
      startConn = startConn%>% left_join(
        gfa$segments %>% select(to = name, toLN = LN), by = "to"
      ) %>% 
        group_by(geneId, from, fromOrient) %>% 
        filter(toLN == max(toLN)) %>% dplyr::slice(1) %>% 
        group_by(geneId, name = from) %>% 
        summarise(sides = n_distinct(fromOrient[toLN >= minBlastLength | second]), 
                  LN1 = max(toLN[fromOrient == "+"], 0),
                  second1 = any(second[fromOrient == "+"], F),
                  LN2 = max(toLN[fromOrient == "-"], 0),
                  second2 = any(second[fromOrient == "-"], F),
                  .groups = "drop"
                  ) 
      
      startConn = bind_rows(
        startConn,
        gfa$segments %>% filter(start > 0, !name %in% startConn$name) %>% 
          select(geneId, name) %>% 
          mutate(sides = 0, LN1 = 0, LN2 = 0, second1 = F, second2 = F)
      )
      
      startConn = startConn %>% 
        left_join(
          gfa$segments %>% select(name, LN, group), by = "name"
        ) %>% 
        #Only keep the longest start within a group
        group_by(geneId, group) %>% 
        filter(LN == max(LN)) %>%
        #Check if start is longest in whole file
        group_by(geneId) %>% 
        mutate(longest = (LN == max(LN))) %>% ungroup()
        
      startConn = startConn %>% 
        mutate(type = case_when(
          (sides == 2 & (second1 & second2)) | (second1 & LN2 >= minBlastLength) |
            (second2 & LN1 >= minBlastLength) | 
            (LN1 >= minBlastLength & LN2 >= minBlastLength) ~ "noFragments",
          sides > 0  & (second1 | second2) ~ "singleWithBranch",
          sides == 1 & (LN1 >= minBlastLength | LN2 >= minBlastLength) ~ "singleLarge",
          TRUE ~ "fragmentsOnly"
        )) 
      
      #! Get startConnInfo before filter 
      myStats$startConn = startConn %>% group_by(geneId) %>% summarise(
        type = case_when(
          any(type == "noFragments") ~ "noFragments",
          all(type == "fragmentsOnly") ~ "fragmentsOnly",
          TRUE ~ "singleSide"
        ), .groups = "drop"
      ) %>% group_by(type) %>% summarise(nBefore = n(), .groups = "drop") %>% 
        mutate(pipelineId = pipelineId)
      #!
      
      #Extract the kmercounts
      kmerCounts = gfa$segments %>% 
        # mutate(keep  = ! group %in% 
        #          (startConn %>% 
        #             filter(type == "fragmentsOnly", LN < 75) %>% 
        #             pull(group) %>% unique())) %>% 
        select(-sequence, startPerc = start) %>%
        mutate(start = startPerc > 0) %>%
        rename(segmentId = name) %>% group_by(geneId) %>%
        filter(any(start)) %>% ungroup()
      
      #Get the ARG list
      myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      argGenes = dbGetQuery(myConn, "SELECT geneId, clusterNr, gene, subtype, nBases FROM ARG") %>%
        mutate(geneId = as.character(geneId))
      dbDisconnect(myConn)
      
      #Detect the most likely genes
      genesDetected = kmerCounts %>%
        # filter(keep) %>% 
        group_by(geneId) %>% mutate(
          nSeg = sum(start),
          fileDepth = sum(KC) / sum(LN),
          maxLN = max(LN[start]),
          LNsum = sum(LN[start]),
          KCsum = sum(KC[start]),
          startPerc = weighted.mean(startPerc[start], LN[start])
        ) %>%
        filter(start, maxLN > 0) %>% 
        filter(LN == max(LN)) %>% dplyr::slice(1) %>% ungroup() %>%
        left_join(argGenes, by = "geneId") %>%
        select(geneId, clusterNr, nBases, gene, subtype, LN, KC, fileDepth,
               nSeg, LNsum, KCsum, startPerc) %>%
        mutate(startDepth = KC / LN) %>% rowwise() %>%
        mutate(
          cover1 = round(min(1, LN / nBases), 4),
          cover2 = round(min(1, LNsum / nBases), 4)) %>% ungroup() %>% 
        # #Change cutoffs?
        filter(cover1 >= 0.1)
      
      #Determine each GFA graph type and add it to the genesDetected
      genesDetected = genesDetected %>% 
        left_join(
          startConn %>% group_by(geneId) %>% 
            summarise(type = ifelse(any(type == "noFragments"),
                                    "noFragments", "fragmented"), 
                      .groups = "drop"), by = "geneId"
        ) %>%
        mutate(pipelineId = pipelineId, runId = runId,
               cover = ifelse(type == "noFragments", cover1, cover2),
               val = cover * startPerc * KC) %>%
        select(pipelineId, runId, everything()) %>% 
        filter(cover >= minCover)
      
      #! Get the geneInfo before filter
      myStats$genesBeforeFilter = list(
        summary = genesDetected %>% summarise(
          n = n(), avgCover = mean(cover), sdCover = sd(cover), minCover = {{minCover}}) %>% 
          mutate(pipelineId = pipelineId), 
        genes = genesDetected$geneId
      )
      #!
      
      #Split the master GFA into fragments and noFragments
      #***************************************************
      
      #When there is at least one noFragments, ignore the rest
      startConn = startConn %>% group_by(geneId) %>% 
        mutate(x = any(type == "noFragments"), y = max(LN)) %>% 
        ungroup() %>% 
        filter((x & LN == y) | !x) %>% select(-x, -y)
      
      #Fragments are short and do not contain any branches
      fragments = gfa$segments %>% 
        filter(geneId %in% genesDetected$geneId,
               group %in% (startConn %>% filter(type == "fragmentsOnly") %>%
                             pull(group) %>% unique())) %>% pull(name)
      
      fragments = gfa_filterSegments(gfa, fragments)
      
      fragments$segments$old = fragments$segments$name
      
      fragments = gfa_mergeSegments(fragments, extraSummaries = list(
        name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_fragUnitig"),
        geneId = function(x) str_extract(x$name[1], "^\\d+"),
        old = function(x) x$name[x$start > 0][1],
        start = function(x) sum(x$start * x$LN) / sum(x$LN)
      ), suffix = "_start")
      
      #Update the startConn with merged segments
      temp = fragments$segments %>% filter(name != old)
      startConn = bind_rows(
        startConn %>% filter(!name %in% temp$old),
        temp %>% transmute(geneId, name, sides = 0, LN1 = 0, second1 = F,
                           LN2 = 0, second2 = F, LN, type = "fragmentsOnly",
                           old) %>% 
          left_join(startConn %>% select(old = name, group, longest), by = "old")
      ) %>% select(-old)
      
      fragments$segments = fragments$segments %>% select(-old) %>% distinct()
      
      #No fragments are all other files (might contain fragments, but ignored)
      noFragments = gfa$segments %>% 
        filter(
          geneId %in% genesDetected$geneId,
          group %in% (startConn %>% filter(type != "fragmentsOnly") %>% 
            pull(group) %>% unique())) %>% pull(name)
      
      noFragments = gfa_filterSegments(gfa, noFragments)
      
      #! Get the updates before fragment filtering
      myStats$removeDuplicates = list(
        pipelineId = pipelineId,
        nFragBefore = n_distinct(fragments$segments$geneId)
      )

      myStats$fragmentsBefore = list(
        segments = fragments$segments %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start > 0), avgSLN = mean(LN[start > 0]), avgSKC = mean(KC[start > 0]),
            .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId)
      )
      
      myStats$noFragmentsBefore = list(
        segments = noFragments$segments %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start > 0), avgSLN = mean(LN[start > 0]), avgSKC = mean(KC[start > 0]),
            .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId),
        links = noFragments$links %>% group_by(geneId) %>% summarise(
          nLinks = n(),.groups = "drop"
        ) %>% mutate(pipelineId = pipelineId)
      )
      #!
      
      #Find identical fragmented files
      #*******************************
      if(nrow(fragments$segments) > 0 & n_distinct(fragments$segments$geneId) > 1){
        
        #check the overlap of identical segments between files
        simMat = fragments$segments %>%
          mutate(val = 1) %>%
          pivot_wider(sequence, names_from = geneId, values_from = LN,
                      values_fill = 0, values_fn = length) %>% select(-sequence) %>%
          as.matrix()
        
        #Build a similarity matrix
        simMat = t(simMat) %*% (simMat > 0)
        simMat = apply(simMat, 2, function(x) x / max(x))
        
        #Extract all identical files and collapse them keeping the
        #one with the highest score
        geneIds = rownames(simMat)
        clusters = apply(simMat, 2, function(x){
          geneIds[x == 1]
        })
        
        #MIf multiple per group, keep only the best
        if(class(clusters) == "list") {
          
          new = data.frame()
          while(length(clusters) > 0){
            new = bind_rows(
              new, genesDetected %>% filter(geneId %in% clusters[[1]]) %>%
                filter(val == max(val)) %>% dplyr::slice(1)
            )
            
            clusters[clusters[[1]]] = NULL
          }
          
          toKeep = new$geneId
          
        } else {
          
          toKeep = clusters
          
        }
        
        fragments = list(
          segments = fragments$segments %>% filter(geneId %in% toKeep),
          links = fragments$links
        )
        
        if(length(toKeep) > 1){
          
          #Temp files for clustering
          myFile = tempfile()
          myFile_RC = tempfile()
          
          #Write fragments to fasta
          x = fragments$segments
          suppressWarnings(fasta_write(x$sequence, myFile,
                                       x$name, type = "n"))
          
          #Add the reverse complement (usearch does not do that when calc_distmx)
          system(sprintf("%s -fastx_revcomp %s -label_suffix _RC -fastaout %s -quiet",
                         settings["usearch"], myFile, myFile_RC))
          invisible(file.append(myFile, myFile_RC))
          
          #Use cluster_fast to reduce number of segments by grouping in identity clusters
          system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.5 -quiet",
                         settings["usearch"], myFile, myFile))
          
          clusters = read.delim(myFile, header = F) %>% 
            mutate(id1 = str_extract(V1, "^\\d+"), id2 = str_extract(V2, "^\\d+")) %>% 
            pivot_longer(c(id1, id2)) %>% 
            mutate(V1 = ifelse(.$V1 > .$V2, .$V2, .$V1),
                   V2 = ifelse(.$V1 > .$V2, .$V1, .$V2)) %>% 
            mutate(id1 = str_extract(V1, "^\\d+"), id2 = str_extract(V2, "^\\d+"),
                   V1 = str_remove(V1, "_RC$"), V2 = str_remove(V2, "_RC$"),
                   V3 = 1 - V3) %>% 
            select(id1, id2, sim = V3, V1, V2) %>% distinct() 
          
          clusters= bind_rows(
            clusters%>% select(id1, id2, sim, V1, V2),
            clusters%>% mutate(id1 = .$id2, id2 = .$id1, sim,
                               V1 = .$V2, V2 = .$V1)
          ) %>% distinct() %>% filter(id1 != id2, sim == 1)
          
          #Only keep clusters where all start segments match (if multiple)
          clusters = clusters%>% 
            filter(str_detect(V1, "_start$"), str_detect(V2, "_start$")) %>% 
            left_join(
              fragments$segments %>% filter(start > 0) %>% group_by(geneId) %>% 
                summarise(nStart = n()),
              by = c("id1"="geneId")
            ) %>% group_by(id1, id2) %>% filter(n_distinct(V1) == nStart[1]) %>% 
            ungroup() %>% select(id1, id2) %>% distinct()
          
          clusters = bind_rows(
            clusters, clusters %>% mutate(id1 = .$id2, id2 = .$id1)
          ) %>% distinct()
          
          #Genes with no clusters are kept per definition
          toKeep = toKeep[!toKeep %in% clusters$id1]
          
          #If two are identical in start, pick the one with overall highest depth 
          depths = fragments$segments %>% group_by(geneId) %>% 
            summarise(depth = sum(KC) / sum(LN), startLN = sum(LN[start > 0]))
          clusters = clusters%>% 
            left_join(depths, by = c("id1" = "geneId")) %>% 
            left_join(depths, by = c("id2" = "geneId")) %>% 
            filter(startLN.x >= startLN.y) %>% 
            arrange(desc(startLN.x), desc(depth.x))
          
          #Remove all genes that are fully contained within other one and have
          #lower depth
          toAdd = c()
          while(!all(unique(clusters$id1) %in% toAdd)){
            myId = unique(clusters$id1)
            myId = myId[!myId %in% toAdd][1]
            
            toRemove = myId
            l = 0
            while(l < length(toRemove)){
              toRemove = clusters %>% filter(id1 %in% toRemove | id2 %in% toRemove) 
              toRemove = unique(c(toRemove$id1, toRemove$id2))
              l = l + 1
            }
            
            toRemove = toRemove[!toRemove %in% myId]
            clusters= clusters%>% filter(!(id1 %in% toRemove | id2 %in% toRemove) |
                                           id1 == myId)
            
            toAdd = c(toAdd, myId)
          }
          
          #Only the remaining fragmented ARG are kept
          toKeep = c(toKeep, unique(clusters$id1))
          
          # #Update the fragments GFA
          # fragments = list(
          #   segments = fragments$segments %>% filter(geneId %in% toKeep),
          #   links = fragments$links 
          # )
          
        }
      } else {
        
        #There are no fragmented files
        toKeep = unique(fragments$segments$geneId)
        
      }
      
      #! Get the updates after frag filtering
      myStats$removeDuplicates$nFragAfter = length(toKeep)
      myStats$removeDuplicates$nBranchedBefore = n_distinct(noFragments$segments$geneId)
      #!

      #Find identical unfragmented files (start segments excluded)
      #**********************************************
      
      if(nrow(noFragments$segments) > 0){
        
        #check the overlap of identical segments between files
        simMat = noFragments$segments %>%
          filter(!str_detect(name, "_start$")) %>%
          mutate(val = 1) %>%
          pivot_wider(sequence, names_from = geneId, values_from = LN,
                      values_fill = 0) %>% select(-sequence) %>%
          as.matrix()
        
        #Build a similarity matrix
        simMat = t(simMat) %*% (simMat > 0)
        simMat = apply(simMat, 2, function(x) x / max(x))
        
        #Extract all identical files and collapse them keeping the
        #one with the highest score
        geneIds = rownames(simMat)
        clusters = apply(simMat, 2, function(x){
          geneIds[x == 1]
        })
        
        
        new = data.frame()
        while(length(clusters) > 0){
          new = bind_rows(
            new, genesDetected %>% filter(geneId %in% clusters[[1]]) %>%
              filter(val == max(val)) %>% dplyr::slice(1)
          )
          
          clusters[clusters[[1]]] = NULL
        }
        
        
        #Examine near identical unfragmented files
        #*****************************************
        simMat = simMat[colnames(simMat) %in% new$geneId,
                        colnames(simMat) %in% new$geneId]
        
        #Get the start segments and immediate neighbours and generate FASTA
        toAlign = noFragments$segments %>%
          filter(geneId %in% new$geneId, str_detect(name, "_start$")) %>%
          group_by(geneId) %>% filter(LN == max(LN)) %>% ungroup() %>%
          mutate(pos = 0)
        
        toAlign = bind_rows(
          toAlign,
          noFragments$links %>% filter(from %in% toAlign$name | to %in% toAlign$name) %>%
            filter(!(str_detect(from, "_start$") & str_detect(to, "_start$"))) %>%
            mutate(name = ifelse(str_detect(from, "_start$"), to, from),
                   pos = ifelse(str_detect(from, "_start$"), 0, 1),
                   pos2 = ifelse(pos == 0, fromOrient, toOrient),
                   pos = case_when(
                     pos == 0 & pos2 == "+" ~ 1,
                     pos == 0 & pos2 == "-" ~ 2,
                     pos == 1 & pos2 == "+" ~ 2,
                     TRUE ~ 1
                   )) %>%
            select(name, pos)  %>%
            left_join(noFragments$segments, by = "name"))
        
        #Find the distance matrix between the seq with usearch
        distMat = map_df(0:2, function(i){
          
          myFile = tempfile()
          myFile_RC = tempfile()
          
          #Write to fasta
          x = toAlign %>% filter(pos == i)
          suppressWarnings(fasta_write(x$sequence, myFile,
                                       x$name, type = "n"))
          
          #Add the reverse complement (usearch does not do that when calc_distmx)
          system(sprintf("%s -fastx_revcomp %s -label_suffix _RC -fastaout %s -quiet",
                         settings["usearch"], myFile, myFile_RC))
          file.append(myFile, myFile_RC)
          
          #Use cluster_fast to reduce number of segments by grouping in identity clusters
          system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.5 -quiet",
                         settings["usearch"], myFile, myFile))
          
          read.delim(myFile, header = F) %>% mutate(pos = i)
          
        })
        
        #Evaluate similarities again
        distMat = distMat %>%
          mutate(id1 = str_extract(V1, "^\\d+"), id2 = str_extract(V2, "^\\d+")) %>%
          pivot_longer(c(id1, id2)) %>%
          mutate(V1 = ifelse(.$V1 > .$V2, .$V2, .$V1),
                 V2 = ifelse(.$V1 > .$V2, .$V1, .$V2)) %>%
          mutate(id1 = str_extract(V1, "^\\d+"), id2 = str_extract(V2, "^\\d+"),
                 V1 = str_remove(V1, "_RC$"), V2 = str_remove(V2, "_RC$"),
                 V3 = 1 - V3) %>%
          select(id1, id2, sim = V3, pos, V1, V2) %>% distinct()
        
        distMat = bind_rows(
          distMat %>% select(id1, id2, sim, pos, V1, V2),
          distMat %>% select(id1 = id2, id2 = id1, sim, pos, V1, V2)
        ) %>% distinct()
        
        
        distMat = distMat %>%
          group_by(V1, V2) %>% filter(sim == max(sim)) %>% ungroup() %>%
          left_join(noFragments$segments %>%
                      select(name, LN1 = LN, KC1 = KC), by = c("V1" = "name")) %>%
          left_join(noFragments$segments %>%
                      select(name, LN2 = LN, KC2 = KC), by = c("V2" = "name")) %>%
          rowwise() %>%
          mutate(ratio = min(LN1, LN2)) %>% group_by(id1, id2) %>%
          mutate(ratio2 = ratio / sum(ratio)) %>% ungroup()
        
        #Get the similarities with scores for each side as well
        distMat = distMat %>%
          group_by(id1, id2) %>% summarise(
            sim = sum(sim * ratio2),
            sim0 = max(sim[pos == 0], 0, na.rm = T),
            sim1 = max(sim[pos == 1], 0, na.rm = T),
            sim2 = max(sim[pos == 2], 0, na.rm = T),
            .groups = "drop"
          ) %>% left_join(
            new %>%
              select(geneId, gene, subtype, val), by = c("id2" = "geneId")
          )
        
        #Filter and group for similarities
        distMat = distMat %>% group_by(id1) %>%
          filter(sim > 0.9, sim0 > 0.9) %>%
          filter(val == max(val)) %>% ungroup() %>%
          select(geneId = id2, gene, subtype) %>% distinct()
        
        # #Update the fragments GFA
        # noFragments = list(
        #   segments = noFragments$segments %>% filter(geneId %in% distMat$geneId),
        #   links = noFragments$links %>% filter(geneId %in% distMat$geneId)
        # ) 
        
        toKeep = c(toKeep, distMat$geneId)
      }
      
      
      #Final filter
      genesDetected = genesDetected %>% filter(geneId %in% toKeep) 
      
      #! Get the genes after filtering
      myStats$removeDuplicates$nBranchedAfter = n_distinct(distMat$geneId)
      myStats$genesAfterFilter = list(
        summary = genesDetected %>% summarise(
          n = n(), avgCover = mean(cover), sdCover = sd(cover), 
          minCover = {{minCover}}, pipelineId = pipelineId), 
        genes = genesDetected$geneId
      )
      #!
     
      #---- Save detected results - but delete old first ----
      myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      
      q = dbSendStatement(myConn, "DELETE FROM detectedARG WHERE pipelineId == ?",
                          params = pipelineId)
      dbClearResult(q)

      q = dbSendStatement(
        myConn,
        paste("INSERT INTO detectedARG ",
              "(pipelineId,runId,geneId,fileDepth,nSeg,LNsum,KCsum,startPerc,LNmax,KCmax,startDepth,cover1,cover2,type)",
              "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"),
              params = unname(as.list(genesDetected %>% 
                  select(pipelineId,runId,geneId,fileDepth,nSeg,LNsum,KCsum,
                         startPerc,LNmax = LN,KCmax = KC,startDepth,cover1,
                         cover2,type))))
      dbClearResult(q)
      dbDisconnect(myConn)

      myDir = sprintf("%sgenesDetected", tempFolder)
      if(dir.exists(myDir)){
        unlink(myDir, recursive = T)
      }
      dir.create(myDir, showWarnings = F)
      write_csv(genesDetected,
                paste0(tempFolder, "genesDetected/genesDetected.csv"))

      #Write the startConn table
      startConn = startConn %>% filter(geneId %in% genesDetected$geneId)
      write_csv(startConn, paste0(tempFolder, "segmentsOfInterest.csv"))
      
      #! Get startConnInfo after filter 
      myStats$startConn = myStats$startConn %>% left_join(
        startConn %>% group_by(geneId) %>% summarise(
        type = case_when(
          any(type == "noFragments") ~ "noFragments",
          all(type == "fragmentsOnly") ~ "fragmentsOnly",
          TRUE ~ "singleSide"
        ), .groups = "drop"
      ) %>% group_by(type) %>% summarise(nAfter = n(), .groups = "drop"),
      by  ="type")
      #!

      #Update GFA files 
      noFragments = list(
        segments = noFragments$segments %>% filter(geneId %in% genesDetected$geneId),
        links = noFragments$links %>% filter(geneId %in% genesDetected$geneId)
      ) 
      
      fragments = list(
        segments = fragments$segments %>% filter(geneId %in% genesDetected$geneId),
        links = fragments$links 
      ) 
      
      #!
      myStats$fragmentsAfter = list(
        segments = fragments$segments %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start > 0), avgSLN = mean(LN[start > 0]), avgSKC = mean(KC[start > 0]),
            .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId)
      )
      
      myStats$noFragmentsAfter = list(
        segments = noFragments$segments %>% 
          group_by(geneId) %>% summarise(
            nSeg = n(), avgLN = mean(LN), sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
            nStart = sum(start > 0), avgSLN = mean(LN[start > 0]), avgSKC = mean(KC[start > 0]),
            .groups = "drop"
          ) %>% mutate(pipelineId = pipelineId),
        links = noFragments$links %>% group_by(geneId) %>% summarise(
          nLinks = n(),.groups = "drop"
        ) %>% mutate(pipelineId = pipelineId)
      )
      #!
      
      #Write files with more than just fragments to separate gfa files 
      #for viewing in Bandage
      for(myGene in unique(noFragments$segments$geneId)){
        myGFA = list()
        myGFA$segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId)
        myGFA$links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId)
        gfa_write(myGFA, sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
      }

      #Write the fragmented ones as one GFA and fasta
      gfa_write(fragments, paste0(tempFolder, "fragmentGFA.gfa"), verbose = 0)

      #Feedback and Logs
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, "Finished detecting ARG"))

    }
  }


  if(maxStep > 2 & nrow(genesDetected) > 0){

    # ---- 3. Simplify the unfragmented or branched GFA files  ---
    #*************************************************************
    if(nrow(logs %>% filter(actionId %in% c(12,14,15))) > 0 & !forceRedo){

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "Skip GFA simplification, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 12,
                                    "Skip GFA simplification, already done"))

      blastSegments = read.csv(sprintf("%sblastSegments.csv", tempFolder))

    } else {

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Simplify GFA files ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 13, "Start simplifying GFA files"))

      dir.create(sprintf("%sgenesDetected/simplifiedGFA", tempFolder), showWarnings = F)
      
      registerDoParallel(cores=maxCPU)

      #Read all GFA files
      #! add .combine = "bind_rows" when removed and chage below too!!
      blastSegments = foreach(myGene = unique(noFragments$segments$geneId)) %dopar% {
                
            fullGFA = gfa_read(sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
            
            #! Get the info before simplification
            tempStats = fullGFA$segments %>% 
              summarise(
                geneId = {{myGene}}, when = "before", 
                nSeg = n(), avgLN = mean(LN), 
                sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
                links = nrow(fullGFA$links), .groups = "drop"
              )
            #!
            
            if(nrow(fullGFA$links) == 0){
              return(data.frame())
            }
            
            segmentsOfInterest = startConn %>% 
              filter(geneId == myGene, type != "fragmentsOnly") %>% 
              pull(name)
            
            myGFA = map(segmentsOfInterest, function(segmentOfInterest){

              #Stay within maxPathDist / maxPathSteps around this segment
              myGFA = gfa_neighbourhood(fullGFA, segmentOfInterest, maxPathDist,
                                        maxSteps = maxPathIter)
              
              #Check if filter yields any results
              if(nrow(myGFA$links) == 0){
                if(verbose > 0){cat("done\n")}
                return(data.frame())
              }
              
              #Remove parallel sequences by picking the shortest
              myGFA = gfa_removeRedundant(myGFA, segmentOfInterest)
              
              myGFA = gfa_mergeSegments(
                myGFA, exclude = segmentOfInterest,
                extraSummaries = list(
                  name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig")
                ))
              
              #Keep pruning the graph to get rid of small appendages
              myGFA = gfa_trimLooseEnds(
                myGFA, trimLength - 1, keepRemoving = T, exclude = segmentOfInterest,
                extraSummaries = list(
                  name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig")
                ))
              
              #Colour the ARG segments green for easier display in Bandage and save results
              myGFA = gfa_annotation(myGFA, myGFA$segments$name[
                str_detect(myGFA$segments$name, "_start$")], color = "green")
              
              return(myGFA)
              
            })
            
            #Merge the GFAs if needed
            if(length(myGFA) > 1){
              fullGFA = list()
              curMaxIndex = 0
              for(i in 1:length(myGFA)){
                
                partGFA = myGFA[[i]]
                unitigs = partGFA$segments$name[
                  str_detect(partGFA$segments$name,"_unitig\\d+$")]
                maxIndex = max(str_extract(unitigs, "\\d+$") %>% as.integer(), 0)
                
                #Update indices if needed
                if(curMaxIndex > 0 & maxIndex > 0){
                  new = paste0(str_extract(unitigs, "\\d+"), "_unitig", 
                               1:length(unitigs) + curMaxIndex)
                  partGFA = gfa_updateNames(partGFA, unitigs, new)
                }
                
                curMaxIndex = curMaxIndex + maxIndex
                
                fullGFA$segments = bind_rows(
                  fullGFA$segments, partGFA$segments %>% mutate(
                    start = replace_na(start, 0), 
                    group = replace_na(group, max(group, na.rm = T))
                  )
                )
                fullGFA$links = bind_rows(fullGFA$links, partGFA$links)
                
              }
            } else {
              
              fullGFA = myGFA[[1]]
              fullGFA$segments = fullGFA$segments %>%  mutate(
                start = replace_na(start, 0), 
                group = replace_na(group, max(group, na.rm = T))
              )
              
            }
            
            gfa_write(fullGFA, sprintf("%s/genesDetected/simplifiedGFA/%s_simplified.gfa",
                                     tempFolder, myGene))
            
            #! Get the info after simplification
            tempStats = bind_rows(
              tempStats, fullGFA$segments %>% 
                summarise(
                  geneId = {{myGene}}, when = "after", 
                  nSeg = n(), avgLN = mean(LN), 
                  sdLN = sd(LN), avgKC = mean(KC), sdKC = sd(KC),
                  links = nrow(fullGFA$links), .groups = "drop"
                )) %>% mutate(pipelineId = pipelineId)
            #!
            
            return(list(gfa = fullGFA$segments %>% 
                          filter((LN > minBlastLength) | (start != "0" & LN >= 75)) %>%
                          mutate(geneId = myGene), tempStats = tempStats))

      }

      #! Get the info after simplification - reduction in cases more than 3 seg
      # myStats$simplification = bind_rows(lapply(blastSegments, "[[", 2)) %>%
      #   group_by(geneId) %>% filter(nSeg[when == "before"] > 3) %>%  summarise(
      #     segReduction = nSeg[when == "after"] / nSeg[when == "before"],
      #     avgLNincrease = avgLN[when == "after"] / avgLN[when == "before"],
      #     sdLNincrease = sdLN[when == "after"] / sdLN[when == "before"],
      #     linkReduction = links[when == "after"] / links[when == "before"],
      #     .groups = "drop"
      #   ) %>% mutate(pipelineId = pipelineId)
      
      myStats$simplification = bind_rows(lapply(blastSegments, "[[", 2)) 
      #!
      
      blastSegments = bind_rows(lapply(blastSegments, "[[", 1))
      
      if(verbose > 0){cat("done\n")}

      #Create final table of fragmented and simplified GFA
      if(nrow(blastSegments) > 0){

        blastSegments = bind_rows(
          blastSegments %>% select(-start, -CL) %>%
            mutate(geneId = as.character(geneId), group = as.numeric(group)),
          fragments$segments %>%
            filter(LN >= minBlastLength | (str_detect(name, "_start$") & LN > 74)) %>%
          distinct()) %>%
          mutate(blastId = name)

      } else {

        blastSegments = fragments$segments %>%
          mutate(blastId = name)

      }

      fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder),
                  blastSegments$blastId, type = "nucleotide")
      
      write.csv(blastSegments, sprintf("%sblastSegments.csv", tempFolder), row.names = F)
      
      #! Blast segments 
      myStats$blastSeg = blastSegments %>% summarise(
        when = "beforeCluster", n = n(), avgLN = mean(LN), sdLN = sd(LN), 
        minBlastLength = {{minBlastLength}},
        .groups = "drop"
      ) 
      #!

      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14,
                                    "Finished simplifying GFA files"))

    }
  } else if (nrow(genesDetected) == 0 & maxStep > 0){
	#Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
						"There were no genes detected. Skip the rest of the pipeline\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 15,
                                    "No genes detected"))

      myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)
      q = dbSendStatement(
        myConn,
        "UPDATE pipeline SET statusCode = 10, statusMessage = 'No ARG detected - pipeline halted', modifiedTimestamp = ?
         WHERE pipelineId = ?",
        params = list(as.character(Sys.time()), pipelineId))
      dbClearResult(q)
      dbDisconnect(myConn)
  }


  if(maxStep > 3 & nrow(genesDetected) > 0){

    # ---- 4. Extract segments for BLAST  ---
    #****************************************
    if(nrow(logs %>% filter(actionId %in% c(16,18))) > 0 & !forceRedo){

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "Skip clustering segments and fasta generation, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 16,
                                    "Skip clustering segments and fasta generation, already done"))

      nFiles = length(list.files(tempFolder, pattern = "blastSegmentsClustered\\d+.fasta"))

    } else {

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "Cluster segments and generate FASTA for BLAST ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 17,
                                    "Start clustering segments and generate FASTA for BLAST"))
      
      #Remove old blast data
      unlink(list.files(tempFolder, ".csv.gz|expand_", full.names = T))

      #Use cluster_fast to reduce number of segments by grouping in identity clusters
      system(sprintf("%s -cluster_fast %s -sort length -query_cov 0.975 -target_cov 0.975 -id %f -uc %s%s -quiet",
                     settings["usearch"],
                     sprintf("%sblastSegments.fasta", tempFolder),
                     clusterIdentidy,
                     sprintf("%sblastSegments.out", tempFolder),
                     ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))

      fastaPaths = read.table(sprintf("%sblastSegments.out", tempFolder))
      fastaPaths = blastSegments %>% filter(blastId %in% (fastaPaths %>% filter(V10 == "*") %>% pull("V9")))
      
      #! Blast segments after cluster
      myStats$blastSeg = bind_rows(myStats$blastSeg, fastaPaths %>% summarise(
        when = "afterCluster", n = n(), avgLN = mean(LN), sdLN = sd(LN), 
        minBlastLength = {{minBlastLength}},
        .groups = "drop"
      )) %>% mutate(pipelineId = pipelineId)
      #!
     
      #FASTA to blast should not contain more than 500,000 nucleotides per file (split if needed)
      x = file.remove(list.files(sprintf("%s",tempFolder),
                             pattern = "blastSegmentsClustered", full.names = T))
      nFiles = ceiling(sum(fastaPaths$LN) / 500000)
      maxPerFile = ceiling(nrow(fastaPaths) / nFiles)

      for(i in 1:nFiles){
        subSet = (1 + maxPerFile*(i-1)):min((i*maxPerFile), nrow(fastaPaths))
        fasta_write(fastaPaths$sequence[subSet],
                    sprintf("%sblastSegmentsClustered%i.fasta", tempFolder, i),
                    fastaPaths$blastId[subSet], type = "n")
      }
      
      #Feedback and Logs
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 18,
                                    "Finished clustering segments and generate FASTA for BLAST"))
    }

  }


  if(maxStep > 4 & nrow(genesDetected) > 0){

    # ---- 5. Finalise preparation ---
    #*********************************

    if(nrow(logs %>% filter(actionId %in% c(16,18))) > 0 & !forceRedo){

      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "No new FASTA files to prepare for BLAST, already done\n\n",
                          "Everything has been successfully run already.\n",
                          " Set forceRedo = TRUE and run again if needed\n\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 19, "No new files to prepare for BLAST"))

    } else {

      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                          "Updating database with new files to BLAST ... ")}

      blastSubmissions = data.frame()
      for(i in 1:nFiles){

        #Prepare the table to insert into blastSubmissions
        blastSubmissions = rbind(
          blastSubmissions,
          list(RID = "",
               timeStamp = as.integer(Sys.time()),
               tempName = tempName,
               fastaFile = sprintf("blastSegmentsClustered%i.fasta", i),
               statusCode = 0,
               statusMessage = "Awaiting submission",
               folder = tempFolder))
      }

      #Update DB
      blastSubmissions$pipelineId = pipelineId
      blastSubmissions$runId = runId
      blastSubmissions = blastSubmissions %>%
        select(pipelineId,runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)

      #Delete old ones fist in case of a redo
      myConn = dbConnect(SQLite(), database, synchronous = NULL)
      sqliteSetBusyHandler(myConn, 30000)

      q = dbSendStatement(myConn, "DELETE FROM blastSubmissions WHERE pipelineId == ?",
                          params = pipelineId)
      dbClearResult(q)

      #Add the new ones
      q = dbSendStatement(
        myConn,
        paste("INSERT INTO blastSubmissions",
              "(pipelineId,runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)",
              "VALUES (?,?,?,?,?,?,?,?,?)"),
        params = unname(as.list(blastSubmissions)))
      dbClearResult(q)

      q = dbSendStatement(
        myConn,
        "UPDATE pipeline SET statusCode = 3, statusMessage = 'Finished blast prep', modifiedTimestamp = ? WHERE pipelineId = ?",
        params = list(as.character(Sys.time()), pipelineId))
      dbClearResult(q)
      dbDisconnect(myConn)

      #Feedback and Logs
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 20,
                                    "Updated database with new files to BLAST"))

    }


  }


  if(maxStep < 5 & nrow(genesDetected) > 0){
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 21,
                                  paste("BLAST prep limited to step", maxStep)))
  } else if(maxStep > 0) {
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 22, "Finished BLAST prep"))
  }

},
finally = {

  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
                      "Cleaning up and saving logs ... ")}

  myConn = dbConnect(SQLite(), database, synchronous = NULL)
  sqliteSetBusyHandler(myConn, 30000)

  #Submit the logs, even in case of error so we know where to resume
  newLogs$runId = runId
  newLogs$tool = "blastPrep.R"
  newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
  q = dbSendStatement(
    myConn,
    "INSERT INTO logs (runId,tool,timeStamp,actionId, actionName) VALUES(?,?,?,?,?)",
    params = unname(as.list(newLogs)))
  dbClearResult(q)
  dbDisconnect(myConn)

  if(file.exists(sprintf("%smasterGFA.db", tempFolder))){
    system(sprintf("%s -f %s", zipMethod, sprintf("%smasterGFA.db", tempFolder)))
  }
  
  #! Wrap up and save 
  myStats$final = data.frame(timeStamp = as.numeric(Sys.time()), 
                             success = 22 %in% newLogs$actionId)
  write_rds(myStats, sprintf("/mnt/meta2amrData/pipelineTest/after1200/stats/%s_stats.rds", pipelineId))
  #!

  #Garbage collection
  x = gc(verbose = FALSE, full = T)

  if(verbose > 0){cat("done\n")}
})
