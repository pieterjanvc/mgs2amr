#!/usr/bin/env Rscript

#**********************
# ---- Blast prep ----
#*********************
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(parallel))

# ---- FUNCTIONS ----
#********************
#Make cut-offs based on steepest tangent in curve (maybe not best if very different RA)
cutOff = function(numbers){
  if(length(numbers) < 2){
    return(numbers)
  }
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique()
}


# ---- Inputs ----
#*****************
args = commandArgs(trailingOnly = TRUE) #arguments specified in the shell file

baseFolder = formatPath(args[[1]], endWithSlash = T)
tempFolder = formatPath(args[[2]], endWithSlash = T)
tempName = args[[3]]
verbose = args[[4]]
runId = as.integer(args[[5]])
pipelineId = as.integer(args[[6]])
keepAllMetacherchantData = as.logical(args[[7]])
maxPathDist = as.integer(args[[8]]) #Distance from ARG to crop the GFA file (reduces blast search)
minBlastLength = as.integer(args[[9]]) #Min segment length to submit to blast
trimLength = as.integer(args[[10]]) #Loose segments smaller than this will be cut from thr GFA
clusterIdentidy  = as.numeric(args[[11]]) #The cluster identity percent used in usearch
forceRedo = as.logical(args[[12]]) #If parts of the code have successfully run before a crash, do not repeat unless forceRedo = T
maxStep= as.integer(args[[13]]) #Which parts of the script to run? If NA all is run

maxStep = ifelse(maxStep == 0, 5, maxStep)

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
myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
logs = dbGetQuery(
  myConn, 
  paste("SELECT * FROM logs WHERE runId IN",
        "(SELECT runId FROM scriptUse WHERE pipelineId == ?)",
        "AND tool = 'blastPrep.R'"), 
  params = pipelineId)
dbDisconnect(myConn)
newLogs = data.frame(timeStamp = as.integer(Sys.time()), 
                     actionId = 1, actionName = "Start BLAST prep")

tryCatch({
  
  if(maxStep > 0){
    
    # ---- 1. Clean up files and folders ----
    #****************************************
    if(nrow(logs %>% filter(actionId %in% c(2, 6))) > 0){
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Skip MetaCherchant cleanup, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                    "Skip MetaCherchant cleanup, already done"))
      
      #Takes long time to load, only do if next step is not completed (not needed afterwards)
      if(nrow(logs %>% filter(actionId %in% c(9,11))) == 0 | forceRedo){
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                            "Load master GFA file for processing ... ")}
        # gfa = gfa_read(gzfile(paste0(tempFolder, "/masterGFA.gfa.gz")))
        system(sprintf("%s -d %s", zipMethod, paste0(tempFolder, "masterGFA.db.gz")))
        myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
        gfa = list(
          segments = dbReadTable(myConn, "segments"),
          links = dbReadTable(myConn, "links")
        )
        dbDisconnect(myConn)
        
        # gfa$segments = read_csv(paste0(tempFolder, "masterGFA_segments.csv.gz")) %>% 
        #   as.data.frame()
        # gfa$links = read_csv(paste0(tempFolder, "masterGFA_links.csv.gz")) %>% 
        #   as.data.frame()
        
        if(verbose > 0){cat("done\n")}
      }
      
    } else { #Not processed yet ...
      
      # ---- Merge all GFA files----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Filter and merge MetaCherchant output ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, "Start filtering and merging MetaCherchant output"))
      
      myFiles = list.files(sprintf("%s", tempFolder), 
                           ".gfa", recursive = T, full.names = T)
      myFiles = myFiles[!str_detect(myFiles, "masterGFA")]

      #This process can be done in parallel so speed things up
      cl <- parallel::makeCluster(detectCores())
      x = clusterEvalQ(cl, {
        library(dplyr)
        library(stringr)
        library(gfaTools)
      })
      clusterExport(cl, varlist = c("tempFolder", "zipMethod"))
      
      #Read all GFA files
      gfa = suppressWarnings(parLapply(cl, myFiles, function(x){
        
        geneId = str_extract(x, "\\d+(?=/graph.gfa)")
        gfa = gfa_fixMetacherchant(x)
        
        #Check if the file is not empty
        if(nrow(gfa$segments) > 0){
          
          gfa$segments$geneId = geneId
          
          if(nrow(gfa$links) > 0){
            gfa$links$geneId = geneId
          }
          
          return(gfa)

        } else {
          return(NULL)
        }
        
      }))
      
      rm(cl)
      
      #Merge all the GFAs that we need for further analysis
      x = sapply(gfa, is.null)
      notUsed = gfa[x] %>% unlist()
      gfa = gfa[!x]
      gfa = list(
        segments = bind_rows(sapply(gfa, "[", 1)),
        links = bind_rows(sapply(gfa, "[", 2))
      )
      rm(x)
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                                    "Finished filtering and merging MetaCherchant output"))
      
      # ---- Update the master GFA file and zip it to save space ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Write master GFA ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Start writing master GFA"))
      
      # gfa_write(gfa, paste0(tempFolder, "masterGFA.gfa"))
      # system(sprintf("%s -f %s", zipMethod, paste0(tempFolder, "masterGFA.gfa")))
      
      myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
      dbWriteTable(myConn, "segments", gfa$segments, overwrite = T)
      dbWriteTable(myConn, "links", gfa$links, overwrite = T)
      dbDisconnect(myConn)
      
      #Remove Metacherchant Data if set
      if(keepAllMetacherchantData == F & (nrow(logs %>% filter(actionId == 4)) == 0)){
        allDirs = list.dirs(tempFolder,recursive = F)
        allDirs = allDirs[!str_detect(allDirs,"metacherchant_logs")]
        system(sprintf("rm -R %s", paste(allDirs, collapse = " ")))
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
      singleSeg = gfa_read(paste0(tempFolder, "fragmentGFA.gfa"))

      
    } else {
      
      # ---- Clean up and merge start segments ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Recover ARG seed sequences ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Start recovering ARG seed sequences"))
      
      #Split the file in groups to merge start segments
      myGroups = gfa$links %>% group_by(geneId) %>% summarise(n = n())
      
      myGroup = 1
      curSum = 0
      myResult = rep(0, nrow(myGroups))
      for(i in 1:nrow(myGroups)){
        
        curSum = curSum + myGroups$n[i]
        
        if(curSum >= 250000){
          curSum = 0
          myGroup = myGroup + 1
        }
        
        myResult[i] = myGroup
      }
      
      myGroups$group = myResult
      
      # myGoups = gfa$segments$geneId
      # steps = c(gfa$segments$geneId[seq(1, length(myGoups), by = 100000)], 
      #           gfa$segments$geneId[nrow(gfa$segments)]) %>% unique()
      # myGoups = unique(gfa$segments$geneId)
      # 
      # if(length(steps) > 2){
      #   myGoups = mapply(function(x, y, z){
      #     z[which(z == x):which(z == y)]
      #   }, x = steps[-length(steps)], y = lead(steps)[-length(steps)], z = list(myGoups))
      #   
      # } else {
      #   myGoups = list(c(myGoups))
      # }
      
      #Run the mergeStartSegments function per group 
      #TODO consider parallel
      # gfa = lapply(myGoups, function(myId){
      #   myGoup = list()
      #   myGoup$segments = gfa$segments %>% filter(geneId %in% myId)
      #   myGoup$links = gfa$links %>% filter(geneId %in% myId)
      #   mergeStartSegments(myGoup, maxGap = 800, maxStep = 20)
      # })
      
      #This process can be done in parallel so speed things up
      cl <- parallel::makeCluster(detectCores())
      x = clusterEvalQ(cl, {
        library(dplyr)
        library(gfaTools)
        library(RSQLite)
      })
      clusterExport(cl, varlist = c("myGroups", "tempFolder"))
      
      #Read all GFA files
      gfa = suppressWarnings(parLapply(cl, unique(myGroups$group), function(myGroup){
        
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
        
        #Add the geneId to the segment names to make sure they are unique 
        gfaGroup$segments = gfaGroup$segments %>% 
          mutate(name = paste0(geneId, "_", name))
        gfaGroup$links = gfaGroup$links %>% 
          mutate(
            from = paste0(geneId, "_", from),
            to = paste0(geneId, "_", to)
          )
        
        #Cut out small appendages from graphs to make joining easier
        gfaGroup = gfa_trimLooseEnds(gfaGroup, 100, keepRemoving = F)
        
        #Merge the start segments
        mergeStartSegments(
          list(
            segments = gfaGroup$segments %>% filter(geneId %in% myGenes) ,
            links = gfaGroup$links %>% filter(geneId %in% myGenes)
          ))
        
        }))
      
      rm(cl)
      
      #Merge the results
      gfa = list(
        segments = bind_rows(sapply(gfa, "[", 1)),
        links = bind_rows(sapply(gfa, "[", 2))
      )
      
     
      #Add the geneId back to the links and unitig names
      gfa$links$geneId = str_extract(gfa$links$from, "^\\d+")

      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, 
                                    "Finished recovering ARG seed sequences"))
      
      
      # ---- Estimate ARG presence by coverage and fragmentation----
      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Detect ARG in the data ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, "Start detecting ARG"))
      
      #Extract the kmercounts
      kmerCounts = gfa$segments %>% select(-sequence, startPerc = start) %>%  
        mutate(start = str_detect(name, "_start$")) %>% 
        rename(segmentId = name)
      
      #Get the ARG list
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
      argGenes = dbGetQuery(myConn, "SELECT geneId, clusterNr, gene, subtype, nBases FROM ARG") %>%
        mutate(geneId = as.character(geneId))

      #Detect the most likely genes
      genesDetected = kmerCounts %>% 
        group_by(geneId) %>% mutate(
          nSeg = sum(start),
          fileDepth = sum(KC) / sum(LN),
          LNsum = sum(LN[start]),
          KCsum = sum(KC[start])
        ) %>% ungroup() %>% 
        filter(start) %>% 
        left_join(argGenes, by = "geneId") %>% 
        group_by(geneId, clusterNr, nBases, gene, subtype, fileDepth, 
                 nSeg, LNsum, KCsum, startPerc) %>% 
        summarise(LNmax = max(LN), KCmax = max(KC), 
                  startPerc = startPerc[LN == LNmax],
                  startDepth = KCmax / LNmax,
                  .groups = 'drop') %>% rowwise() %>% 
        # filter(gene == "blaTEM") %>%
        mutate(
          cover1 = round(min(1, LNmax / nBases), 4),
          cover2 = round(min(1, LNsum / nBases), 4),
          val = cover1 * startPerc * KCmax) %>% 
        group_by(gene) %>% 
        mutate(simToBest = between(fileDepth, 
                                   fileDepth[val == max(val)][1] * 0.9,
                                   fileDepth[val == max(val)][1] * 1.1)) %>% 
        group_by(gene, simToBest) %>% arrange(desc(val)) %>% slice(1) %>% 
        ungroup() %>% 
        filter(cover1 >= 0.25, startDepth / fileDepth >= 0.25) %>% 
        select(-val, -simToBest)
        # filter(adjCov == max(adjCov)) %>% 
        # group_by(gene) %>% filter(
        #   startPerc == max(startPerc) |
        #     cover1 == max(cover1)) %>% 
      # group_by(clusterNr, gene, subtype) %>% 
      #   filter(cover1 == max(cover1), startDepth == max(startDepth)) %>% 
      #   ungroup() %>% 
      #   mutate(pipelineId = pipelineId, runId = runId) %>% 
      #   group_by(gene, subtype) %>% filter(cover1 == max(cover1)) %>% slice(1) %>% 
      #   group_by(clusterNr) %>% filter(cover1 == max(cover1)) %>% slice(1) %>% 
      #   ungroup() %>% arrange(desc(cover1))

      #Only keep genes that are minimum 90% covered (or use cut-off when higher) 
      # genesDetected = genesDetected %>% 
      #   filter(cover1 >= cutOff(genesDetected$cover1))
      
      # ---- Get all fragmented GFA (only start seg or one connection) ----
      #********************************************************************
      #Get all start segments
      singleSeg = gfa$segments %>% 
        filter(
          geneId %in% genesDetected$geneId,
          str_detect(name, "_start$"))
      
      mySegments = singleSeg$name
      
      #Find all the segments that connect to a start segment
      singleSeg = gfa$links %>% 
        filter(
          geneId %in% genesDetected$geneId,
          from %in% mySegments | 
            to %in% mySegments) 
      
      #Save start segments that do not connect to anything (isolated)
      mySegments = mySegments[!mySegments %in% unique(c(singleSeg$from, singleSeg$to))]
      
      singleSeg = c(singleSeg$from, singleSeg$to) %>% unique()
      singleSeg = singleSeg[!str_detect(singleSeg, "_start$")]
      
      singleSeg = gfa$links %>%
        filter(
          geneId %in% genesDetected$geneId,
          from %in% singleSeg |
            to %in% singleSeg) %>%
        mutate(start = (str_detect(from, "_start") |
                          str_detect(to, "_start")))
      
      singleSeg = singleSeg %>% group_by(geneId) %>%
        filter(!any(
          from[start] %in% c(from[!start], to[!start])) &
            !any(to[start] %in% c(from[!start], to[!start])))
      
      if(nrow(singleSeg) > 0){
        singleSeg = singleSeg %>% select(from, to, geneId)
        singleSeg[(nrow(singleSeg)+1):(nrow(singleSeg)*2),] = data.frame(
          from  = singleSeg$to,
          to = singleSeg$from,
          geneId = singleSeg$geneId
        )
      }
       
      
      singleSeg = singleSeg %>% 
        # distinct() %>% 
        filter(str_detect(from, "_start$")) %>% 
        group_by(from) %>% summarise(n = n()) %>% 
        filter(n < 3) %>% 
        pull(from)
      
      
      
      # # singleSeg = c(singleSeg$from, singleSeg$to) %>% unique()
      # # 
      # # #Save start segments that do not connect to anything (isolated)
      # # mySegments = mySegments[!mySegments %in% singleSeg]
      # # 
      # #Get all segments that only connect to a start segment (i.e. are end segments)
      # singleSeg = singleSeg[!str_detect(singleSeg, "_start$")]
      # test = c(singleSeg$from, singleSeg$to) %>% unique()
      # test = test[!str_detect(test, "_start$")]
      # 
      # test = gfa$links %>%
      #   filter(
      #     geneId %in% genesDetected$geneId,
      #     from %in% test |
      #       to %in% test) %>%
      #   mutate(start = (str_detect(from, "_start") |
      #                     str_detect(to, "_start")))
      # 
      # singleSeg = test %>% group_by(geneId) %>%
      #   filter(!any(
      #     from[start] %in% c(from[!start], to[!start])) &
      #       !any(to[start] %in% c(from[!start], to[!start])))
      # 
      
      #GFA should ONLY have single fragments - not combi ...
      
      singleSeg = unique(c(mySegments, singleSeg))
      
      #Add graphs that only consist of a single start segment
      onlyStartSeg = gfa$segments$geneId[!gfa$segments$geneId %in% gfa$links$geneId]
      onlyStartSeg = onlyStartSeg[onlyStartSeg %in% genesDetected$geneId]
      
      #Build a GFA that contains all these short GFAs (will be individual islands)
      singleSeg = list(segments = gfa$segments %>% 
                         filter(name %in% singleSeg | geneId %in% onlyStartSeg), 
                       links = gfa$links %>% filter(
                         from %in% singleSeg | to %in% singleSeg))
      
      #Consider fragment only if small piece of gene 
      myFilter = gfa$segments %>% select(name, LN, geneId) %>% 
        filter(geneId %in% singleSeg$segments$geneId, 
               name %in% c(singleSeg$links$from, singleSeg$links$to)) %>% 
        group_by(geneId) %>% summarise(LN = sum(LN), n = n()) %>% 
        filter(!(n < 4 & LN > 500))
      
      singleSeg = list(segments = gfa$segments %>% 
                         filter(name %in% singleSeg | geneId %in% onlyStartSeg), 
                       links = gfa$links %>% filter(
                         from %in% singleSeg | to %in% singleSeg))
      
      #Check fragment type
      fragType = gfa$segments %>% select(name, geneId, LN) %>% 
        filter(geneId %in% singleSeg$segments$geneId) %>% 
        filter(str_detect(name, "_start$"))
      
     
      
      if(nrow(fragType) > 0){
        fragType = fragType %>% left_join(singleSeg$segments %>% select(-sequence, -LN),
                                          by = c("name", "geneId")) %>% 
          group_by(geneId) %>% 
          summarise(type = case_when(
            sum(is.na(KC)) == 0 ~ "fragmentsOnly",
            is.na(sum(KC[LN == max(LN)])) ~ "longestFragment",
            TRUE ~ "longestNoFragment"
          ), .groups = "drop")
      } else {
        fragType = data.frame(geneId = "", type = NA)
      }
      
      
      
      #Add to this new info to detected genes table
      genesDetected = genesDetected %>% left_join(fragType, by = "geneId") %>% 
        mutate(type = replace_na(type, "noFragments")) %>% 
        #Only keep decently covered genes (assumed higher abundance), or fragments (low abundance)
        filter((cover1 > 0.1 & type != "fragmentsOnly") | type == "fragmentsOnly") %>% 
        mutate(pipelineId = pipelineId, runId = runId) %>% 
        select(pipelineId, runId, everything())
      
      #---- Save detected results - but delete old first ----
      q = dbSendStatement(myConn, "DELETE FROM detectedARG WHERE pipelineId == ?", 
                          params = pipelineId)
      dbClearResult(q)
      
      q = dbSendStatement(
        myConn, 
        paste("INSERT INTO detectedARG ",
              "(pipelineId,runId,geneId,fileDepth,nSeg,LNsum,KCsum,startPerc,LNmax,KCmax,startDepth,cover1,cover2,type)",
              "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"), 
              params = unname(as.list(
                genesDetected %>% select(-gene, -subtype, -clusterNr, -nBases))))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      dir.create(sprintf("%sgenesDetected", tempFolder), showWarnings = F)
      write_csv(genesDetected, 
                paste0(tempFolder, "genesDetected/genesDetected.csv"))
      
      #Write the unfragmented gfa files as separate files for viewing in Bandage
      for(myGene in genesDetected$geneId[! genesDetected$type %in% 
                                         c("fragmentsOnly", "longestFragment")]){
        myGFA = list()
        myGFA$segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId)
        myGFA$links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId)
        gfa_write(myGFA, sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
      }
      
      #Write the fragmented ones as a one GFA and fasta
      gfa_write(singleSeg, paste0(tempFolder, "fragmentGFA.gfa"), verbose = 0)
      # gfa_writeUnitigs(singleSeg, paste0(tempFolder, "blastSegFragment.fasta"))
      
      #Feedback and Logs
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, "Finished detecting ARG"))
      
    }
  }
  
  
  if(maxStep > 2 & nrow(genesDetected) > 0){
    
    # ---- 3. Simplify the unfragmented GFA files  ---
    #*************************************************
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
      if(verbose > 0){cat("\n")}
      
      blastSegments = map_df(
        genesDetected$geneId[! genesDetected$type %in% 
                               c("fragmentsOnly", "longestFragment")], function(myGene){

        if(verbose > 0){
          cat(sprintf(" gene %i/%i ... ", 
                      which(myGene == genesDetected$geneId[! genesDetected$type %in% 
                                                             c("fragmentsOnly", "longestFragment")]),
                      length(genesDetected$geneId[! genesDetected$type %in% 
                                                    c("fragmentsOnly", "longestFragment")])))}
        
        gfa = gfa_read(sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
        
        #Check if filter yields any results
        if(nrow(gfa$links) == 0){
          return(data.frame())
        }
        
        #Get largest start segment
        segmentOfInterest = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
          filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% slice(1) %>% pull(name)
        
        #Stay within maxPathDist around this segment
        gfa = gfa_neighbourhood(gfa, segmentOfInterest, maxPathDist, noLoops = T)
        
        #Check if filter yields any results
        if(nrow(gfa$links) == 0){
          if(verbose > 0){cat("done\n")}
          return(data.frame())
        }
        
        #Get all other start segments too
        allStart = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
          pull(name)
        
        #Remove parallel sequences by picking the shortest
        gfa = gfa_trimLooseEnds(gfa, trimLength, keepRemoving = F, exclude = segmentOfInterest)
        gfa = gfa_removeRedundant(gfa, segmentOfInterest, maxLN = trimLength)
        gfa = gfa_mergeSegments(
          gfa, exclude = segmentOfInterest,
          extraSummaries = list(
            name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig")
          ))
        
        #Keep pruning the graph to get rid of small appendages
        gfa = gfa_trimLooseEnds(
          gfa, trimLength - 1, keepRemoving = T, exclude = segmentOfInterest,
          extraSummaries = list(
            name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig")
          ))
        
        #Colour the ARG segments green for easier display in Bandage and save results
        gfa = gfa_annotation(gfa, gfa$segments$name[
          str_detect(gfa$segments$name, "_start$")], color = "green")
        gfa_write(gfa, sprintf("%s/genesDetected/simplifiedGFA/%s_simplified.gfa",
                               tempFolder, myGene))
        
        if(verbose > 0){cat("done\n")}
        
        gfa$segments %>% filter(LN > minBlastLength) %>% 
          mutate(geneId = myGene)
        
      }) 
      
      id = genesDetected$geneId[genesDetected$type == "fragmentsOnly"]
      fragmentsOnly = list() 
      fragmentsOnly$segments = gfa$segments %>% filter(geneId %in% id)
      fragmentsOnly$links = gfa$links %>% filter(geneId %in% id)

      fragmentsOnly = gfa_mergeSegments(fragmentsOnly, extraSummaries = list(
        name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig"),
        geneId = function(x) str_extract(x$name[1], "^\\d+")
      ))
      
      gfa_write(fragmentsOnly, sprintf("%sfragmentsOnly.gfa", tempFolder))
      
      #Create final table of fragmented and simplified GFA
      if(nrow(blastSegments) > 0){

        blastSegments = bind_rows(
          blastSegments %>% select(-start, -CL), 
          fragmentsOnly$segments %>% select(-start) %>% 
            filter(LN >= 100)) %>% distinct() %>% 
          mutate(blastId = name) 
        
      } else {
        
        blastSegments = fragmentsOnly$segments %>% 
          mutate(blastId = name)
        
      }
      
      fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder), 
                  blastSegments$blastId, type = "nucleotide")
      write.csv(blastSegments, sprintf("%sblastSegments.csv", tempFolder), row.names = F)
      
      #Feedback and Logs
      if(verbose > 0){cat(" finished\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14, 
                                    "Finished simplifying GFA files"))
      
    }
  } else if (nrow(genesDetected) == 0){
	#Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),
						"There were no genes detected. Skip the rest of the pipeline\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 15, 
                                    "No genes detected"))
      
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
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
      
      #Use cluster_fast to reduce number of segments by grouping in identity clusters
      system(sprintf("%s -cluster_fast %s -sort size -query_cov 0.75 -id %f -uc %s%s",
                     settings["usearch"],
                     sprintf("%sblastSegments.fasta", tempFolder),
                     clusterIdentidy,
                     sprintf("%sblastSegments.out", tempFolder),
                     ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))
      
      fastaPaths = read.table(sprintf("%sblastSegments.out", tempFolder))
      fastaPaths = blastSegments %>% filter(blastId %in% (fastaPaths %>% filter(V10 == "*") %>% pull("V9")))
      
      #FASTA to blast should not contain more than 250,000 nucleotides per file (split if needed)
      nFiles = ceiling(sum(fastaPaths$LN) / 250000)
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
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
      
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
  } else {
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 22, "Finished BLAST prep"))
  }
  
}, 
finally = {
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
  #Submit the logs, even in case of error so we know where to resume
  newLogs$runId = runId
  newLogs$tool = "blastPrep.R"
  newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
  q = dbSendStatement(
    myConn, 
    "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
    params = unname(as.list(newLogs)))
  dbClearResult(q)
  dbDisconnect(myConn)  
  
  if(file.exists(sprintf("%smasterGFA.db", tempFolder))){
    system(sprintf("%s -f %s", zipMethod, sprintf("%smasterGFA.db", tempFolder)))
  }
  
})