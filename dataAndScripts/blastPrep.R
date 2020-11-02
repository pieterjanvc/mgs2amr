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
keepAllMetacherchantData = as.logical(args[[6]])
maxPathDist = as.integer(args[[7]]) #Distance from ARG to crop the GFA file (reduces blast search)
minBlastLength = as.integer(args[[8]]) #Min segment length to submit to blast
trimLength = as.integer(args[[9]]) #Loose segments smaller than this will be cut from thr GFA
clusterIdentidy  = as.numeric(args[[10]]) #The cluster identity percent used in usearch
forceRedo = as.logical(args[[11]]) #If parts of the code have successfully run before a crash, do not repeat unless forceRedo = T
# prevRunId = as.integer(args[[12]]) #If present, use info from prev run to skip part that were completed

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
prevRunId = ifelse(file.exists(paste0(tempFolder,"runId")),
                   readLines(paste0(tempFolder,"runId"), n = 1) %>% as.integer(), 0)

#Check the log file to see if there was a previous run of the code
myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
logs = dbGetQuery(myConn, "SELECT * FROM logs WHERE runId = ? AND tool = 'blastPrep.R'", 
                  params = prevRunId)
dbDisconnect(myConn)
newLogs = data.frame(timeStamp = as.integer(Sys.time()), 
                     actionId = 1, actionName = "Start BLAST prep")

tryCatch({
  # ---- Clean up files and folders ----
  #*************************************
  if(nrow(logs %>% filter(actionId %in% c(2, 4))) > 0){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                        "Skip MetaCherchant cleanup, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                  "Skip MetaCherchant cleanup, already done"))
  
    #Takes long time to load, only do if next step is not completed (not needed afterwards)
    if(nrow(logs %>% filter(actionId %in% c(5,7))) == 0){
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                          "Load master GFA file for processing ... ")}
      gfa = gfa_read(gzfile(paste0(tempFolder, "masterGFA.gfa.gz"), 
                            "masterGFA.gfa"))
      if(verbose > 0){"done"}
    }
    
  } else {
    
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
    
    #Update the master GFA file and zip it to save space
    gfa_write(gfa, paste0(tempFolder, "masterGFA.gfa"))
    system(sprintf("%s %s", paste0(tempFolder, "masterGFA.gfa"), zipMethod))
    
    #Remove Metacherchant Data if set
    if(keepAllMetacherchantData == F & (nrow(logs %>% filter(actionId == 4)) == 0)){
      allDirs = list.dirs(tempFolder,recursive = F)
      allDirs[!str_detect(allDirs,"metacherchant_logs")]
      system(sprintf("rm -R %s", paste(allDirs, collapse = " ")))
    }
    
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, "Finished merging MetaCherchant output"))
  
  }
  
  
  
  # ---- Detect important ARG ----
  #*******************************
  if(nrow(logs %>% filter(actionId %in% c(5, 7))) > 0 & !forceRedo){
    
    #Update the runId for the detectedARG in the database
    myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
    q = dbSendStatement(myConn, "UPDATE detectedARG SET runId = ? WHERE runId = ?", 
                        params = list(runId, prevRunId))
    dbClearResult(q)
    dbDisconnect(myConn)
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Skip ARG detection, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Skip ARG detection, already done"))
    
    genesDetected = read.csv(paste0(tempFolder, "genesDetected/genesDetected.csv"))
    
  } else {
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Detect ARG in the data ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, "Start detecting ARG"))
    
    #Extract the kmercounts
    kmerCounts = gfa$segments %>% 
      mutate(start = str_detect(name, "_start$")) %>% 
      rename(segmentId = name)
    
    #Get the ARG list
    myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
    argGenes = dbGetQuery(myConn, "SELECT geneId, clusterNr, nBases FROM ARG") %>%
      mutate(geneId = as.character(geneId))
    
    #Detect the most likely genes
    genesDetected = kmerCounts %>% filter(start) %>% 
      left_join(argGenes, by = c("geneId" = "geneId")) %>% 
      group_by(geneId, clusterNr, nBases) %>% 
      summarise(segmentLength = sum(LN), kmerCount = sum(KC), n = n(), .groups = 'drop') %>% rowwise() %>% 
      mutate(coverage = round(min(1, segmentLength / nBases)), 4) %>% 
      group_by(clusterNr) %>% 
      filter(kmerCount == max(kmerCount)) %>% ungroup() %>% 
      mutate(runId = runId, geneId = as.integer(geneId)) %>% 
      select(runId, geneId, segmentLength, kmerCount, n, coverage)
    
    
    #Only keep genes that are minimum 90% covered (or use cut-off when higher) 
    genesDetected = genesDetected %>% 
      filter(coverage >= max(0.9, cutOff(genesDetected$covered)))
    
    #Save results 
    q = dbSendStatement(myConn, "INSERT INTO detectedARG VALUES(?,?,?,?,?,?)", 
                        params = unname(as.list(genesDetected)))
    dbClearResult(q)
    dbDisconnect(myConn)
    
    dir.create(sprintf("%sgenesDetected", tempFolder), showWarnings = F)
    write.csv(genesDetected %>% select(-runId), 
              paste0(tempFolder, "genesDetected/genesDetected.csv"), row.names = F)
    
    #Write the detected gfa files as separate files for view in Bandage
    for(myGene in genesDetected$geneId){
      myGFA = list()
      myGFA$segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId)
      myGFA$links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId)
      gfa_write(myGFA, sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
    }
    
    #Feedback and Logs
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Finished detecting ARG"))
    
  }
  
  
  
  # ---- Simplify the GFA files  ---
  #*************************************
  if(nrow(logs %>% filter(actionId %in% c(8, 10))) > 0 & !forceRedo){
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Skip GFA simplification, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, "Skip GFA simplification, already done"))
    
    blastSegments = read.csv(sprintf("%sblastSegments.csv", tempFolder))
    
  } else {
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Simplify GFA files ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Start simplifying GFA files"))
    
    dir.create(sprintf("%sgenesDetected/simplifiedGFA", tempFolder), showWarnings = F)
    if(verbose > 0){cat("\n")}
    
    # #Don't repeat any simplifications that have been done in the past if doing it again
    # alreadyDone = list.files(sprintf("%sgenesDetected/simplifiedGFA", tempFolder),
    #                          pattern = "_simplified.gfa") %>% str_extract("^\\d+")
    
    blastSegments = map_df(genesDetected$geneId, function(myGene){
      
      # if(verbose > 0){cat(sprintf(" %s (%s) ...", genesDetected[genesDetected$geneId == myGene, "ncbi"],
      #                             genesDetected[genesDetected$geneId == myGene, "gene"]))}
      if(verbose > 0){cat(sprintf(" gene %i/%i ... ", which(myGene == genesDetected$geneId), 
                                  length(genesDetected$geneId)))}
      
      gfa = gfa_fixMetacherchant(sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
      
      #Start from the longest ARG segment (might need to be updated in future)
      startSegment = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
        filter(LN == max(LN)) %>% pull(name)
      
      #Check if filter yields any results
      if(nrow(gfa$links) == 0){
        return(data.frame())
      }
      
      #Detect al paths to the ARG and only keep those segments
      allPaths = gfa_pathsToSegment(gfa, startSegment[1], maxDistance = maxPathDist, 
                                    verbose = F, segmentListOnly = T)
      allPaths = sapply(allPaths, "[[", "segmentOrder") %>% unlist %>% unique()
      gfa = gfa_filterSegments(gfa, allPaths, action = "keep")
      
      #Check if filter yields any results
      if(nrow(gfa$links) == 0){
        return(data.frame())
      }
      
      #Trim loose small segments and merge segments in series
      for(i in 1:3){ # TODO create better loop
        gfa = gfa_trimLooseEnds(gfa, trimLength, verbose = F)
      }
      gfa = gfa_mergeSegments(gfa, gfa$segments$name[str_detect(gfa$segments$name, "_start$")])
      
      #Colour the ARG segments green for easier display in Bandage and save results
      gfa = gfa_annotation(gfa, gfa$segments$name[str_detect(gfa$segments$name, "_start$")], color = "green")
      gfa_write(gfa, sprintf("%s/genesDetected/simplifiedGFA/%s_simplified.gfa", tempFolder, myGene))
      
      if(verbose > 0){cat("done\n")}
      
      gfa$segments %>% filter(LN > minBlastLength, !str_detect(name, "_start$")) %>% 
        mutate(geneId = myGene)
      
    }) %>% mutate(blastId = paste0(geneId, "_", name))
    
    #Write a FASTA file with all segments that should be submitted to BLAST
    fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder),
                blastSegments$blastId, type = "n")
    
    write.csv(blastSegments, sprintf("%sblastSegments.csv", tempFolder), row.names = F)
    
    #Feedback and Logs
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, "Finished simplifying GFA files"))
                        
  }
  
  
  # ---- Extract segments for BLAST  ---
  #*************************************
  if(nrow(logs %>% filter(actionId %in% c(11, 13))) > 0 & !forceRedo){
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                        "Skip clustering segments and fasta generation, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, "Skip clustering segments and fasta generation, already done"))
    
    nFiles = length(list.files(tempFolder, pattern = "blastSegmentsClustered\\d+.fasta"))
    
  } else {
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), "Cluster segments and generate FASTA for BLAST ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 12, "Start clustering segments and generate FASTA for BLAST"))
    
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
      fasta_write(fastaPaths$sequence[subSet], sprintf("%sblastSegmentsClustered%i.fasta", tempFolder, i),
                  fastaPaths$blastId[subSet], type = "n")
    }
    
    #Feedback and Logs
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 13, "Finished clustering segments and generate FASTA for BLAST"))
  }
  
  
  # ---- Finalise preparation ---
  #******************************
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
  if(nrow(logs %>% filter(actionId %in% c(14, 15))) > 0 & !forceRedo){
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                        "No new FASTA files to prepare for BLAST, already done\n\n",
                        "Everything has been successfully run already.\n",
                        " Set forceRedo = TRUE and run again if needed\n\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14, "No new files to prepare for BLAST"))
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"),"Updating database with new files to BLAST ... ")}
    
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
    blastSubmissions$runId = runId
    blastSubmissions = blastSubmissions %>% 
      select(runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)
    q = dbSendStatement(myConn, paste("INSERT INTO blastSubmissions",
                                      "(runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)",
                                      "VALUES (?,?,?,?,?,?,?,?)"), 
                        params = unname(as.list(blastSubmissions)))
    dbClearResult(q)
    
    #Feedback and Logs
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 15, "Updated database with new files to BLAST"))
    
  }
  
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 16, "Finished BLAST prep"))
}, 
finally = {
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
  #Submit the logs, even in case of error so we know where to resume
  newLogs$runId = runId
  newLogs$tool = "blastPrep.R"
  newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
  q = dbSendStatement(myConn, "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
                      params = unname(as.list(newLogs)))
  dbClearResult(q)
  
  dbDisconnect(myConn)
  
  #Add the runId to the temp folder
  write(runId, paste0(tempFolder, "runId"))
  
})