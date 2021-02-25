#!/usr/bin/env Rscript

#**********************
# ---- Blast prep ----
#*********************
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RSQLite))

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
    if(nrow(logs %>% filter(actionId %in% c(2, 4))) > 0){
      
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Skip MetaCherchant cleanup, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                    "Skip MetaCherchant cleanup, already done"))
      
      #Takes long time to load, only do if next step is not completed (not needed afterwards)
      if(nrow(logs %>% filter(actionId %in% c(5,7))) == 0 | forceRedo){
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                            "Load master GFA file for processing ... ")}
        gfa = gfa_read(gzfile(paste0(tempFolder, "/masterGFA.gfa.gz")))
        if(verbose > 0){"done"}
      }
      
    } else { #Not processed yet ...
      
      # ---- Merge all GFA files----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Merge MetaCherchant output ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, "Start merging MetaCherchant output"))
      
      if(nrow(logs %>% filter(actionId == 4)) == 0){
        #Merge all metacherchant output GFAs
        system(sprintf("cat $(find %s -name '*.gfa') > %smasterGFA.gfa", tempFolder, tempFolder))
      }
      
      #Read the master GFA file as gfa object (fix MetaCherchant error!)
      filePath = paste0(tempFolder, "masterGFA.gfa")
      gfa = gfa_fixMetacherchant(filePath)
      
      #Read the raw file
      myFile = str_split(readLines(filePath), "\t")
      
      #Get all IDs for which the file is not empty
      geneId = system(paste("find ", tempFolder, " -name '*.gfa' | xargs wc -l"), intern = T)
      
      
      geneId = str_match(geneId, "(\\d+).*/(\\d+)/graph\\.gfa$")
      geneId = data.frame(geneId = geneId[,3], count = as.integer(geneId[,2])) %>% 
        na.omit() 
      
      #Add the gene ID to the segments
      gfa$segments$geneId = 
        rep(geneId$geneId, geneId$count)[which(sapply(myFile, "[[", 1) == "S")]
      
      #Add the gene ID to the links
      gfa$links$geneId = 
        rep(geneId$geneId, geneId$count)[which(sapply(myFile, "[[", 1) == "L")]
      
      rm(myFile)
      
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                                    "Finished merging MetaCherchant output"))
      
      # ---- Update the master GFA file and zip it to save space ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Write master GFA to zip ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Start writing master GFA to zip"))
      
      gfa_write(gfa, paste0(tempFolder, "masterGFA.gfa"))
      system(sprintf("%s %s", zipMethod, paste0(tempFolder, "masterGFA.gfa")))
      
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
    if(nrow(logs %>% filter(actionId %in% c(5, 7))) > 0 | forceRedo){
      
      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Skip ARG detection, already done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Skip ARG detection, already done"))
      
      genesDetected = read.csv(paste0(tempFolder, "genesDetected/genesDetected.csv"))
      if(nrow(logs %>% filter(actionId %in% c(8,10))) == 0 | forceRedo){
        singleSeg = gfa_read(paste0(tempFolder, "fragmentGFA.gfa"))
      }
      
    } else {
      
      # ---- Clean up and merge start segments ----
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Recover ARG seed sequences ... ")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Start recovering ARG seed sequences"))
      
      #Add the geneId to the segment names to make sure they are unique 
      gfa$segments = gfa$segments %>% 
        mutate(name = paste0(geneId, "_", name))
      gfa$links = gfa$links %>% 
        mutate(
          from = paste0(geneId, "_", from),
          to = paste0(geneId, "_", to)
        )
      
      #Perform the start merging algorithm on all data at once
      gfa = mergeStartSegments(gfa, maxGap = 500, maxStep = 25)
      
      #Add the geneId back to the links
      gfa$links = gfa$links %>% left_join(
        gfa$segments %>% select(name, geneId),
        by = c("from" = "name")
      )
      
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
        mutate(
          cover1 = round(min(1, LNmax / nBases), 4),
          cover2 = round(min(1, LNsum / nBases), 4)) %>% 
        group_by(gene) %>% filter(startPerc == max(startPerc)) %>% 
      group_by(clusterNr, gene, subtype) %>% 
        filter(cover1 == max(cover1), startDepth == max(startDepth)) %>% 
        ungroup() %>% 
        mutate(pipelineId = pipelineId, runId = runId) %>% 
        group_by(gene, subtype) %>% filter(cover1 == max(cover1)) %>% slice(1) %>% 
        group_by(clusterNr) %>% filter(cover1 == max(cover1)) %>% slice(1) %>% 
        ungroup() %>% arrange(desc(cover1))

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
      singleSeg = c(singleSeg$from, singleSeg$to) %>% unique()
      
      #Save start segments that do not connect to anything (isolated)
      mySegments = mySegments[!mySegments %in% singleSeg]
      
      #Get all segments that only connect to a start segment (i.e. are end segments)
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
      
      #GFA should ONLY have single fragments - not combi ...
      
      singleSeg = unique(c(mySegments, singleSeg$from, singleSeg$to))
      
      #Add graphs that only consist of a single start segment
      onlyStartSeg = gfa$segments$geneId[!gfa$segments$geneId %in% gfa$links$geneId]
      
      #Build a GFA that contains all these short GFAs (will be individual islands)
      singleSeg = list(segments = gfa$segments %>% 
                         filter(name %in% singleSeg | geneId %in% onlyStartSeg), 
                       links = gfa$links %>% filter(
                         from %in% singleSeg | to %in% singleSeg))
      
      #Check fragment type
      fragType = gfa$segments %>% select(name, geneId, LN) %>% 
        filter(geneId %in% singleSeg$segments$geneId) %>% 
        filter(str_detect(name, "_start$"))
      
      fragType = fragType %>% left_join(singleSeg$segments %>% select(-sequence, -LN),
                         by = c("name", "geneId")) %>% 
        group_by(geneId) %>% 
        summarise(type = case_when(
          sum(is.na(KC)) == 0 ~ "fragmentsOnly",
          is.na(sum(KC[LN == max(LN)])) ~ "longestNoFragment",
          TRUE ~ "longestFragment"
        ))
      
      #Add to this new info to detected genes table
      genesDetected = genesDetected %>% left_join(fragType, by = "geneId") %>% 
        mutate(type = replace_na(type, "noFragments")) %>% 
        #Only keep decently covered genes (assumed higher abundance), or fragments (low abundance)
        filter((cover1 > 0.5 & type != "fragmentsOnly") | type == "fragmentsOnly") %>% 
        select(pipelineId, runId, everything())
      
      # genesDetected$fragmented = genesDetected$geneId %in% singleSeg$segments$geneId
      
      #---- Save detected results - but delete old first ----
      q = dbSendStatement(myConn, "DELETE FROM detectedARG WHERE pipelineId == ?", 
                          params = pipelineId)
      dbClearResult(q)
      
      q = dbSendStatement(
        myConn, 
        paste("INSERT INTO detectedARG ",
              "(pipelineId,runId,geneId,fileDepth,nSeg,LNsum,KCsum,startPerc,LNmax,KCmax,startDepth,cover1,cover2,type)",
              "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"), 
              params = unname(as.list(
                genesDetected %>% select(-gene, -subtype, -clusterNr, -nBases))))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      dir.create(sprintf("%sgenesDetected", tempFolder), showWarnings = F)
      write.csv(genesDetected, 
                paste0(tempFolder, "genesDetected/genesDetected.csv"), row.names = F)
      
      #Write the unfragmented gfa files as separate files for viewing in Bandage
      for(myGene in genesDetected$geneId[genesDetected$type != "fragmentsOnly"]){
        myGFA = list()
        myGFA$segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId)
        myGFA$links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId)
        gfa_write(myGFA, sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
      }
      
      #Write the fragmented ones as a one GFA and fasta
      gfa_write(singleSeg, paste0(tempFolder, "fragmentGFA.gfa"))
      # gfa_writeUnitigs(singleSeg, paste0(tempFolder, "blastSegFragment.fasta"))
      
      #Feedback and Logs
      if(verbose > 0){cat("done\n")}
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, "Finished detecting ARG"))
      
    }
  }
  
  
  if(maxStep > 2 & nrow(genesDetected) > 0){
    
    # ---- 3. Simplify the unfragmented GFA files  ---
    #*************************************************
    if(nrow(logs %>% filter(actionId %in% c(8, 10))) > 0 & !forceRedo){
      
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
      
      blastSegments = map_df(genesDetected$geneId[genesDetected$type != "fragmentsOnly"], function(myGene){
        
        if(verbose > 0){cat(sprintf(" gene %i/%i ... ", 
                                    which(myGene == genesDetected$geneId[genesDetected$type != "fragmentsOnly"]), 
                                    length(genesDetected$geneId[genesDetected$type != "fragmentsOnly"])))}
        
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
          return(data.frame())
        }
        
        #Get all other start segments too
        allStart = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
          pull(name)
        
        #Keep pruning the graph to get rid of small appendages
        gfa = gfa_trimLooseEnds(gfa, trimLength, keepRemoving = T, exclude = segmentOfInterest)
        
        #Remove parallel sequences by picking the shortest
        gfa = gfa_removeRedundant(gfa, segmentOfInterest, maxLN = trimLength)
        

        #Colour the ARG segments green for easier display in Bandage and save results
        gfa = gfa_annotation(gfa, gfa$segments$name[
          str_detect(gfa$segments$name, "_start$")], color = "green")
        gfa_write(gfa, sprintf("%s/genesDetected/simplifiedGFA/%s_simplified.gfa", 
                               tempFolder, myGene))
        
        if(verbose > 0){cat("done\n")}
        
        gfa$segments %>% filter(LN > minBlastLength, !str_detect(name, "_start$")) %>% 
          mutate(geneId = myGene)
        
      }) 
      
      #Create final table of fragmented and simplified GFA
      if(nrow(blastSegments) > 0){

        blastSegments = bind_rows(
          blastSegments %>% select(-start), 
          singleSeg$segments %>% select(-start) %>% filter(LN >= 100)) 
        
        #Write a FASTA file with all segments that should be submitted to BLAST
        fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder),
                    blastSegments$blastId, type = "n")
        
      } else {
        
        blastSegments = singleSeg$segments
        
        # #Create an empty file of simplified blast segments
        # system(sprintf("rm %s; touch %s", sprintf("%sblastSegments.fasta", tempFolder),
        #                sprintf("%sblastSegments.fasta", tempFolder)))
      }
      
      blastSegments = blastSegments %>% mutate(blastId = paste0(geneId, '_', name))
      
      # #Add the fragmented segments
      # system(sprintf("cat %sblastSegFragment.fasta %sblastSegments.fasta > %sblastSegments.fasta ", 
      #        tempFolder, tempFolder, tempFolder))
      
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
    if(nrow(logs %>% filter(actionId %in% c(11, 13))) > 0 & !forceRedo){
      
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
      system(sprintf("%s -cluster_fast %s -sort size -id %f -uc %s%s",
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
    
    if(nrow(logs %>% filter(actionId %in% c(14, 15))) > 0 & !forceRedo){
      
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
      
      #Insert into DB
      blastSubmissions$pipelineId = pipelineId
      blastSubmissions$runId = runId
      blastSubmissions = blastSubmissions %>% 
        select(pipelineId,runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)
      
      myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
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
})