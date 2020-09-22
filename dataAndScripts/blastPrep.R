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
usearch = args[[4]]
verbose = args[[5]]
runId = as.integer(args[[6]])
keepAllMetacherchantData = as.logical(args[[7]])
maxPathDist = as.integer(args[[8]]) #Distance from ARG to crop the GFA file (reduces blast search)
minBlastLength = as.integer(args[[9]]) #Min segment length to submit to blast
trimLength = as.integer(args[[10]]) #Loose segments smaller than this will be cut from thr GFA
clusterIdentidy  = as.numeric(args[[11]]) #The cluster identity percent used in usearch
forceRedo = as.logical(args[[12]]) #If parts of the code have successfully run before a crash, do not repeat unless forceRedo = T

# baseFolder = "/data/aplab/ARG_PJ/aim2/meta2amr/"
# tempFolder = "/scratch/van9wf/pipelineTemp/testFile_1599247967_1599248674/"
# tempName = "testFile_1599247967_1599248674"
# usearch =  "/usr/local/usearch/10.0.240/bin/usearch"
# verbose = 1
# keepAllMetacherchantData = T
# maxPathDist = 5000 #Distance from ARG to crop the GFA file (reduces blast search)
# minBlastLength = 250 #Min segment length to submit to blast
# trimLength = 100 #Loose segments smaller than this will be cut from thr GFA
# clusterIdentidy  = 0.95 #The cluster identity percent used in usearch
# forceRedo = F #If parts of the code have successfully run before a crash, do not repeat unless forceRedo = T
tempFolder = formatPath(paste0(tempFolder, tempName), endWithSlash = T)
logPath = sprintf("%s%s_log.csv", tempFolder,tempName)
prevRunId = ifelse(file.exists(paste0(tempFolder,"runId")),
                   readLines(paste0(tempFolder,"runId"), n = 1) %>% as.integer(), 0)

#Check the log file to see if there was a previous run of the code
myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
logs = dbGetQuery(myConn, "SELECT * FROM logs WHERE runId = ? AND tool = 'blastPrep.R'", params = prevRunId)
dbDisconnect(myConn)
newLogs = data.frame(timeStamp = as.integer(Sys.time()), actionId = 1, actionName = "Start BLAST prep")

tryCatch({
  # ---- Clean up files and folders ----
  #*************************************
  if(nrow(logs %>% filter(actionId == 4)) > 0){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Skip MetaCherchant cleanup: Can only be done once (read previous GFA file)\n")}
    
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, "Skip MetaCherchant cleanup: Can only be done once"))
  
    #Takes long time to load, only do if next step is not completed (not needed afterwards)
    if(nrow(logs %>% filter(actionId == 7)) == 0){
      gfa = gfa_read(gzfile(paste0(tempFolder, "masterGFA.gfa.gz"), "masterGFA.gfa"))
    }
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Merge MetaCherchant output ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, "Start merging MetaCherchant output"))
  
    if(nrow(logs %>% filter(actionId == 4)) == 0){
      #Merge all metacherchant output GFAs
      system(sprintf("cat $(find %s -name '*.gfa') > %smasterGFA.gfa", tempFolder, tempFolder))
    }
    
    #Read the master GFA file as gfa object
    filePath = paste0(tempFolder, "masterGFA.gfa")
    gfa = gfa_read(filePath)
    
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
    system(sprintf("gzip %s", paste0(tempFolder, "masterGFA.gfa")))
    
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
  if(nrow(logs %>% filter(actionId == 7)) > 0 & !forceRedo){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Skip ARG detection, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, "Finished merging MetaCherchant output"))
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Skip ARG detection, already done"))
    
    genesDetected = read.csv(paste0(tempFolder, "genesDetected/genesDetected.csv"))
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Detect ARG in the data ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, "Start detecting ARG"))
    
    #Extract the kmercounts
    kmerCounts = gfa$segments %>% 
      mutate(start = str_detect(name, "_start$")) %>% 
      rename(segmentId = name)
    
    #Get the ARG list
    argGenes = read.csv(paste0(baseFolder, "dataAndScripts/argTable.csv"), 
                        colClasses = list(geneId = "character"))
    
    #Detect the most likely genes
    genesDetected = kmerCounts %>% filter(start) %>% 
      left_join(argGenes %>% select(geneId, clusterNr, gene, subtype, nBases, name), by = c("geneId" = "geneId")) %>% 
      group_by(geneId, clusterNr, gene, subtype, nBases, name) %>% 
      summarise(length = sum(LN), kmerCount = sum(KC), n = n(), .groups = 'drop') %>% rowwise() %>% 
      mutate(covered = min(1, length / nBases)) %>% 
      extract(name, into = c("ncbi", "descr"), regex = ">([^\\s]+)\\s+(.*)") %>% 
      group_by(clusterNr) %>% 
      filter(kmerCount == max(kmerCount))
    
    
    #Only keep genes that are minimum 90% covered (or use cut-off when higher) 
    genesDetected = genesDetected %>% filter(covered >= max(0.9, cutOff(genesDetected$covered)))
    
    #Save results 
    dir.create(sprintf("%sgenesDetected", tempFolder), showWarnings = F)
    write.csv(genesDetected, paste0(tempFolder, "genesDetected/genesDetected.csv"), row.names = F)
    
    #Write the detected gfa files as separate files for view in Bandage
    for(myGene in genesDetected$geneId){
      myGFA = list()
      myGFA$segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId)
      myGFA$links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId)
      gfa_write(myGFA, sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
    }
    
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Finished detecting ARG"))
    
  }
  
  
  
  # ---- Simplify the GFA files  ---
  #*************************************
  if(nrow(logs %>% filter(actionId == 10)) > 0 & !forceRedo){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Skip GFA simplification, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, "Skip GFA simplification, already done"))
    
    blastSegments = read.csv(sprintf("%sblastSegments.csv", tempFolder))
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Simplify GFA files ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Start simplifying GFA files"))
    
    dir.create(sprintf("%sgenesDetected/simplifiedGFA", tempFolder), showWarnings = F)
    if(verbose > 1){cat("\n")}
    
    blastSegments = map_df(genesDetected$geneId, function(myGene){
      
      if(verbose > 1){cat(sprintf(" %s (%s) ...", genesDetected[genesDetected$geneId == myGene, "ncbi"],
                                  genesDetected[genesDetected$geneId == myGene, "gene"]))}
      
      gfa = gfa_read(sprintf("%sgenesDetected/%s.gfa", tempFolder, myGene))
      
      #Start from the longest ARG segment (might need to be updated in future)
      startSegment = gfa$segments %>% filter(str_detect(name, "_start$")) %>%
        filter(LN == max(LN)) %>% pull(name)
      
      #Limit the GFA to neighbourhood around the ARG set by maxPathDist
      gfa = gfa_neighbourhood(gfa, startSegment[1], maxPathDist)
      
      #Detect al paths to the ARG and only keep those segments
      allPaths = gfa_pathsToSegment(gfa, startSegment[1], returnList = T)
      allPaths = sapply(allPaths, "[[", "segmentOrder") %>% unlist %>% unique()
      gfa = gfa_filterSegments(gfa, allPaths, action = "keep")
      
      #Trim loose small segments and merge segments in series
      for(i in 1:3){ # TODO create better loop
        gfa = gfa_trimLooseEnds(gfa, trimLength, verbose = F)
      }
      gfa = gfa_mergeSegments(gfa, gfa$segments$name[str_detect(gfa$segments$name, "_start$")])
      
      #Colour the ARG segments green for easier display in Bandage and save results
      gfa = gfa_annotation(gfa, gfa$segments$name[str_detect(gfa$segments$name, "_start$")], color = "green")
      gfa_write(gfa, sprintf("%s/genesDetected/simplifiedGFA/%s_simplified.gfa", tempFolder, myGene))
      
      if(verbose > 1){cat("done\n")}
      
      gfa$segments %>% filter(LN > minBlastLength, !str_detect(name, "_start$")) %>% 
        mutate(geneId = myGene)
      
    }) %>% mutate(blastId = paste0(geneId, "_", name))
    
    #Write a FASTA file with all segments that should be submitted to BLAST
    fasta_write(blastSegments$sequence, sprintf("%sblastSegments.fasta", tempFolder),
                blastSegments$blastId, type = "n")
    
    write.csv(blastSegments, sprintf("%sblastSegments.csv", tempFolder), row.names = F)
    
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, "Finished simplifying GFA files"))
                        
  }
  
  
  # ---- Extract segments for BLAST  ---
  #*************************************
  if(nrow(logs %>% filter(actionId == 13)) > 0 & !forceRedo){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), 
                        "Skip clustering segments and fasta generation, already done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, "Skip clustering segments and fasta generation, already done"))
    
    nFiles = length(list.files(tempFolder, pattern = "blastSegmentsClustered\\d+.fasta"))
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), "Cluster segments and generate FASTA for BLAST ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 12, "Start clustering segments and generate FASTA for BLAST"))
    
    #Use cluster_fast to reduce number of segments by grouping in identity clusters
    system(sprintf("%s -cluster_fast %s -sort size -id %f -uc %s%s",
                   usearch,
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
    
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 13, "Finished clustering segments and generate FASTA for BLAST"))
  }
  
  
  # ---- Submit to BLAST  ---
  #**************************
  if(nrow(logs %>% filter(actionId == 16)) > 0 & !forceRedo){
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "), 
                        "Skip submitting FASTA files to BLAST, already done\n\nEverything has been successfully run already.",
                        "Set forceRedo = TRUE and run again if needed\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14, "Skip submitting FASTA files to BLAST, already done"))
    
  } else {
    
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S   "),"Submitting FASTA files to BLAST ... ")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 15, "Start submitting FASTA files to BLAST"))
    
    blastSubmissions = data.frame()
    for(i in 1:nFiles){
      
      Sys.sleep(3) #Prevent too fast sequential submissions (api will block those)
      
      #Submit to the blast API with limits to bacterial genomes (taxid2)
      # RID = blast_submit(sprintf("%sblastSegmentsClustered%i.fasta", tempFolder, i),
      #                    program = "blastn", megablast = T, database = "nt", expect = 1e-25, word_size = 64, max_num_seq = 50,
      #                    entrez_query = "txid2 [ORGN]", verbose = 1)
      RID = "FAKEID"
      
      #Save the RID info to a central file (following up the submission is not done in this process as it can take a while)
      blastSubmissions = rbind(
        blastSubmissions,
        list(RID = RID,
             timeStamp = as.integer(Sys.time()), 
             tempName = tempName, 
             fastaFile = sprintf("blastSegmentsClustered%i.fasta", i), 
             statusCode = ifelse(is.na(RID), 1, 2), 
             statusMessage = ifelse(is.na(RID), "Failed submission", "Successful submission"), 
             folder = tempFolder))
    }
    
    if(verbose > 0){cat("done\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 16, "Finished submitting FASTA files to BLAST"))
  }
  
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 17, "Finished BLAST prep"))
}, 
finally = {
  
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
  #Submit the info on BALST submissions (if no fail)
  if(exists("blastSubmissions")){
    blastSubmissions$runId = runId
    blastSubmissions = blastSubmissions %>% 
      select(runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)
    q = dbSendStatement(myConn, paste("INSERT INTO blastSubmissions",
                                      "(runId,RID,timeStamp,tempName,fastaFile,statusCode,statusMessage,folder)",
                                      "VALUES (?,?,?,?,?,?,?,?)"), 
                        params = unname(as.list(blastSubmissions)))
    dbClearResult(q)
  }
  
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