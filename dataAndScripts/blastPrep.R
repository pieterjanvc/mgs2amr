#!/usr/bin/env Rscript

#**********************
# ---- Blast prep ----
#*********************
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(parallel))


# ---- Inputs ----
#*****************
args = commandArgs(trailingOnly = TRUE) #arguments specified in the shell file

baseFolder = formatPath(args[[1]], endWithSlash = T)
tempFolder = formatPath(args[[2]], endWithSlash = T)
tempName = args[[3]]
verbose = abs(as.integer(args[[4]]))
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

maxStartGap = 800
 
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
      # if(nrow(logs %>% filter(actionId %in% c(9,11))) == 0 | forceRedo){
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

        if(verbose > 0){cat("done\n")}
      # }
      
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
      
      myConn = dbConnect(SQLite(), sprintf("%smasterGFA.db", tempFolder))
      dbWriteTable(myConn, "segments", gfa$segments, overwrite = T)
      dbWriteTable(myConn, "links", gfa$links, overwrite = T)
      dbDisconnect(myConn)
      
      #Remove Metacherchant Data if set
      if(keepAllMetacherchantData == F & (nrow(logs %>% filter(actionId == 4)) == 0)){
        # allDirs = list.dirs(tempFolder,recursive = F)
        # allDirs = allDirs[!str_detect(allDirs,"metacherchant_logs")]
        # system(sprintf("rm -R '%s'", paste(allDirs, collapse = " ")))
        system(sprintf("ls -d %s/* | grep -P \"/\\d+$\" | xargs rm -R", tempFolder))
        # system(sprintf('find %s -name "lcl*" -exec rm -r {} \\; > /dev/null 2>&1', tempFolder))
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
      
      #Split the file in groups to merge start segments (multiple can be done at once)
      myGroups = gfa$links %>% group_by(geneId) %>% summarise(n = n())
      
      myGroup = 1
      curSum = 0
      myResult = rep(0, nrow(myGroups))
      for(i in 1:nrow(myGroups)){
        
        curSum = curSum + myGroups$n[i]
        
        if(curSum >= 500000){
          curSum = 0
          myGroup = myGroup + 1
        }
        
        myResult[i] = myGroup
      }
      
      myGroups$group = myResult
      
      #This process can be done in parallel to speed things up
      cl <- parallel::makeCluster(min(detectCores(), myGroup))
      x = clusterEvalQ(cl, {
        library(dplyr)
        library(gfaTools)
        library(RSQLite)
        library(stringr)
      })
      clusterExport(cl, varlist = c("myGroups", "tempFolder", "maxStartGap"))
      
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
        
        #Add the geneId to the segment/link names
        gfaGroup$segments = gfaGroup$segments %>% 
          mutate(name = paste0(geneId, "_", name))
        gfaGroup$links = gfaGroup$links %>% 
          mutate(
            from = paste0(geneId, "_", from),
            to = paste0(geneId, "_", to)
          )
        
        #Cut out small appendages from graphs to make joining easier
        gfaGroup = gfa_trimLooseEnds(gfaGroup, 100, keepRemoving = F)
        
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
        
        }))
      stopCluster(cl)
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
        ) %>% 
        filter(start) %>% filter(LN == max(LN)) %>% ungroup() %>% 
        left_join(argGenes, by = "geneId") %>% 
        select(geneId, clusterNr, nBases, gene, subtype, LN, KC, fileDepth, 
               nSeg, LNsum, KCsum, startPerc) %>% 
        mutate(startDepth = KC / LN) %>% rowwise() %>% 
        # filter(gene == "blaTEM") %>%
        mutate(
          cover1 = round(min(1, LN / nBases), 4),
          cover2 = round(min(1, LNsum / nBases), 4),
          val = cover1 * startPerc * KC) %>% ungroup() %>% 
        filter((fileDepth >= 10 & cover1 > 0.75) | fileDepth < 10)
      
      
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

      #GFA should ONLY have single fragments - not combi ...
      
      singleSeg = c(mySegments, singleSeg) %>% 
        str_extract("^\\d+") %>% unique()
      
      #If seed is well covered, it's not a fragment but linear assembly
      singleSeg = genesDetected %>% 
        filter(geneId %in% singleSeg, cover1 < 0.2) %>% pull(geneId)
      
      # singleSeg$segments %>% 
      #   filter(!str_detect(name, "_start$")) %>% 
      #   group_by(geneId) %>% 
      #   filter(n() < 3 & sum(LN) > 750) %>% pull(geneId) %>% unique()
      # 
      # #Add graphs that only consist of a single start segment
      # onlyStartSeg = gfa$segments$geneId[!gfa$segments$geneId %in% gfa$links$geneId]
      # onlyStartSeg = onlyStartSeg[onlyStartSeg %in% genesDetected$geneId]
      # 
      # #Build a GFA that contains all these short GFAs (will be individual islands)
      # singleSeg = list(segments = gfa$segments %>% 
      #                    filter(name %in% singleSeg | geneId %in% onlyStartSeg), 
      #                  links = gfa$links %>% filter(
      #                    from %in% singleSeg | to %in% singleSeg))
      # 
      # #Consider fragment only if small piece of gene 
      # myFilter = gfa$segments %>% select(name, LN, geneId) %>% 
      #   filter(geneId %in% singleSeg$segments$geneId, 
      #          name %in% c(singleSeg$links$from, singleSeg$links$to)) %>% 
      #   group_by(geneId) %>% summarise(LN = sum(LN), n = n()) %>% 
      #   filter(!(n < 4 & LN > 500))
      
      singleSeg = list(segments = gfa$segments %>% filter(geneId %in% singleSeg), 
                       links = gfa$links %>% filter(geneId %in% singleSeg))
      
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
        # #Only keep decently covered genes (assumed higher abundance), or fragments (low abundance)
        # filter((cover1 > 0.1 & type != "fragmentsOnly") | type == "fragmentsOnly") %>% 
        mutate(pipelineId = pipelineId, runId = runId) %>% 
        select(pipelineId, runId, everything())
      
      
      
      #Find identical files (start segments excluded)
      #**********************************************
      
      #check the overlap of identical segments between files
      simMat = gfa$segments %>% 
        filter(!str_detect(name, "_start$"),
               geneId %in% genesDetected$geneId) %>% 
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
      test = apply(simMat, 2, function(x){
        geneIds[x == 1]
      })
      
      
      new = data.frame()
      while(length(test) > 0){
        new = bind_rows(
          new, genesDetected %>% filter(geneId %in% test[[1]]) %>% 
            filter(val == max(val)) %>% slice(1)
        )
        
        test[test[[1]]] = NULL
      }
      
      
      #Examine near identical files
      #****************************
      simMat = simMat[colnames(simMat) %in% new$geneId, 
                      colnames(simMat) %in% new$geneId]
      
      #Get the start segments and immediate neighbours and generate FASTA
      toAlign = gfa$segments %>% 
        filter(geneId %in% new$geneId, str_detect(name, "_start$")) %>% 
        group_by(geneId) %>% filter(LN == max(LN)) %>% ungroup() %>% 
        mutate(pos = 0)
      toAlign = bind_rows(
        toAlign,
        gfa$links %>% filter(from %in% toAlign$name | to %in% toAlign$name) %>% 
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
          left_join(gfa$segments, by = "name"))
      
      #Find the distance matrix between the seq with usearch
      distMat = map_df(0:2, function(i){
        
        myFile = tempfile()
        myFile_RC = tempfile()
        
        #Write to fasta
        x = toAlign %>% filter(pos == i)
        fasta_write(x$sequence, myFile, 
                    x$name, type = "n")
        
        #Add the reverse complement (usearch does not do that when calc_distmx)
        system(sprintf("%s -fastx_revcomp %s -label_suffix _RC -fastaout %s",
                       settings["usearch"], myFile, myFile_RC))
        file.append(myFile, myFile_RC)
        
        #Use cluster_fast to reduce number of segments by grouping in identity clusters
        system(sprintf("%s -calc_distmx %s -tabbedout %s -termdist 0.5",
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
        left_join(gfa$segments %>% 
                    select(name, LN1 = LN, KC1 = KC), by = c("V1" = "name")) %>% 
        left_join(gfa$segments %>% 
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
      
      #Finally only keep the relevant genes
      genesDetected = genesDetected %>% filter(geneId %in% distMat$geneId)
      
      
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
      
      myDir = sprintf("%sgenesDetected", tempFolder)
      if(dir.exists(myDir)){
        unlink(myDir, recursive = T)
      }
      dir.create(myDir, showWarnings = F)
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

      #Process in parallel
      cl <- parallel::makeCluster(detectCores())
      x = clusterEvalQ(cl, {
        library(dplyr)
        library(stringr)
        library(gfaTools)
        library(igraph)
      })
      clusterExport(cl, 
                    varlist = c("genesDetected", "tempFolder", "verbose",
                                "maxPathDist", "trimLength", "minBlastLength"))
      
      #Read all GFA files
      blastSegments = suppressWarnings(parLapply(cl,
        genesDetected$geneId[! genesDetected$type %in% 
                               c("fragmentsOnly", "longestFragment")], function(myGene){

            gfa = gfa_read(sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
            
            #Check if filter yields any results
            if(nrow(gfa$links) == 0){
              return(data.frame())
            }
            
            #Remove disconnected pieces that do not contain a start segment
            myGraph = graph_from_data_frame(data.frame(
              from = gfa$links$from,
              to = gfa$links$to
            ), directed = F)
            
            myGraph = components(myGraph)$membership
            myGraph = data.frame(
              name = names(myGraph),
              group = myGraph
            ) %>%
              group_by(group) %>%
              filter(any(str_detect(name, "_start$"))) %>% 
              ungroup() %>% distinct()
            
            gfa = gfa_filterSegments(gfa, myGraph$name[myGraph$name %in% gfa$segments$name])
            
            #Get largest start segment
            segmentOfInterest = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
              filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% slice(1) %>% pull(name)
            
            #Stay within maxPathDist / maxPathSteps around this segment
            gfa = gfa_neighbourhood(gfa, segmentOfInterest, maxPathDist)
            
            #Check if filter yields any results
            if(nrow(gfa$links) == 0){
              if(verbose > 0){cat("done\n")}
              return(data.frame())
            }
            
            #Get all other start segments too
            allStart = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
              pull(name)
            
            #Remove parallel sequences by picking the shortest
            gfa = gfa_removeRedundant(gfa, segmentOfInterest)
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
        
        return(gfa$segments %>% filter((LN > minBlastLength) | 
                                         (start == "1" & LN >= 75)) %>% 
          mutate(geneId = myGene))
        
      })) 
      stopCluster(cl)
      rm(cl)
      blastSegments = bind_rows(blastSegments)
      if(verbose > 0){cat("done\n")}
      
      #Extract fragmented genes from master GFA
      id = genesDetected$geneId[genesDetected$type == "fragmentsOnly"]
      fragmentsOnly = list() 
      fragmentsOnly$segments = gfa$segments %>% filter(geneId %in% id)
      fragmentsOnly$links = gfa$links %>% filter(geneId %in% id)

      fragmentsOnly = gfa_mergeSegments(fragmentsOnly, extraSummaries = list(
        name = function(x) paste0(str_extract(x$name[1], "^\\d+"), "_unitig"),
        geneId = function(x) str_extract(x$name[1], "^\\d+")
      ))
      
      gfa_write(fragmentsOnly, sprintf("%sfragmentsOnly.gfa", tempFolder), verbose = 0)
      
      #Create final table of fragmented and simplified GFA
      if(nrow(blastSegments) > 0){

        blastSegments = bind_rows(
          blastSegments %>% select(-start, -CL) %>% 
            mutate(geneId = as.character(geneId)), 
          fragmentsOnly$segments %>% 
            filter(LN >= 100 | (str_detect(name, "_start$") & LN > 74)) %>% 
          distinct()) %>% 
          mutate(blastId = name)
        
      } else {
        
        blastSegments = fragmentsOnly$segments %>% 
          mutate(blastId = name)
        
      }
      
      fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder), 
                  blastSegments$blastId, type = "nucleotide")
      write.csv(blastSegments, sprintf("%sblastSegments.csv", tempFolder), row.names = F)
      
      #Feedback and Logs
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
      system(sprintf("%s -cluster_fast %s -sort length -query_cov 0.975 -target_cov 0.975 -id %f -uc %s%s",
                     settings["usearch"],
                     sprintf("%sblastSegments.fasta", tempFolder),
                     clusterIdentidy,
                     sprintf("%sblastSegments.out", tempFolder),
                     ifelse(verbose < 2, " >/dev/null 2>&1", " 2>&1")))
      
      fastaPaths = read.table(sprintf("%sblastSegments.out", tempFolder))
      fastaPaths = blastSegments %>% filter(blastId %in% (fastaPaths %>% filter(V10 == "*") %>% pull("V9")))
      
      #FASTA to blast should not contain more than 250,000 nucleotides per file (split if needed)
      x = file.remove(list.files(sprintf("%s",tempFolder), 
                             pattern = "blastSegmentsClustered", full.names = T))
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
  
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                      "Cleaning up and saving logs ... ")}
  
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
  
  #Garbage collection
  x = gc(verbose = FALSE, full = T)
  
  if(verbose > 0){cat("done\n")}
})
