#!/usr/bin/env Rscript

#******************************
# ---- Local BLASTn search ----
#******************************
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))

#Make sure the blast bin folder is in the R env 
#Sys.getenv("PATH")
#If not, add export PATH=/opt/ncbi-blast-2.13.0+/bin:$PATH to PATH in 
#Renviron file at /opt/R/4.0.2/lib/R/etc/Renviron

# blastn = "/opt/ncbi-blast-2.13.0+/bin/blastn"
# blastDB = "/mnt/meta2amrData/ncbi/blastdb/"

# ---- Inputs ----
#*****************

args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file
baseFolder = formatPath(args[[1]], endWithSlash = T)
database = args[[2]]
runId = as.integer(args[[3]])
verbose = abs(as.integer(args[[4]]))
blastn = args[[5]]
blastDB = args[[6]]
pipelineId = str_trim(unlist(strsplit(args[[7]], ",")))
forceRedo = as.logical(args[[8]])

maxCPU = parallel::detectCores() ## EDIT LATER 

Sys.setenv(BLASTDB = blastDB)

#Check if pigz is available

zipMethod = ifelse(system("command -v pigz", ignore.stdout = T) != 1, 
                   "pigz", "gzip")

#Set these general blast args
blastArgs = list(
  db = "nt",
  evalue = "1e-20",
  word_size = 64,
  max_target_seqs = 250,
  max_hsps = 3,
  taxidlist = sprintf("%sdataAndScripts/%s", baseFolder, "bact.txids"),
  outfmt = "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart send bitscore score length pident nident qcovs qcovhsp"
)


# ---- FUNCTIONS ----
#********************
blast_readOutput = function(file, outfmt, separate = T, includeIssues = F, verbose = 1){
  
  #Load the output csv (zipped)
  blastOut = read.table(file, sep = "\t", quote = "", comment.char = "")
  colnames(blastOut) = strsplit(outfmt, " ")[[1]][-1]
  
  #split multiple matches
  blastOut = blastOut %>% 
    mutate(x = str_count(staxids, ";"), y = str_count(sallacc, ";"), 
           z = x == y | x == 0, plasmid = str_detect(salltitles, "plasmid|Plasmid")) 
  
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


# ---- MAIN CODE ----
#********************
# myQuery = "SELECT * FROM blastSubmissions WHERE"
# 
# if(length(pipelineId) != 0 & !forceRedo){
#   myQuery = paste(myQuery, "statusCode in (0,1,3,4) AND")
# } else {
#   myQuery = paste(myQuery, "statusCode in (0,1,3,4)")
# }
# 
# #Limit the search for specific pipelineIds if set, else do all
# myQuery = paste(myQuery, ifelse(length(pipelineId) != 0,
# 	sprintf("runId in (SELECT runId FROM scriptUse WHERE pipelineId IN ('%s'))", 
# 	        paste(pipelineId, collapse = "','")),
# 	""))

myQuery = sprintf("SELECT * FROM blastSubmissions WHERE %s %s %s", 
        ifelse(forceRedo & length(pipelineId) != 0, "", "statusCode in (0,1,3,4)"), 
        ifelse(!forceRedo & length(pipelineId) != 0, "AND", ""),
        ifelse(length(pipelineId) != 0, 
               sprintf("pipelineId IN ('%s')",
                       paste(pipelineId, collapse = "','")), ""
               ))

#Get all blast submissions that need further follow-up
myConn = dbConnect(SQLite(), database, synchronous = NULL)
sqliteSetBusyHandler(myConn, 30000)
submTable = dbGetQuery(myConn, myQuery) %>% arrange(timeStamp)
dbDisconnect(myConn) 

#Status codes
statusCodes = data.frame(
  status = 0:13,
  message = c("awaiting submission","successful remote submission", "failed remote submission", "searching remotely", 
              "remote search completed", "timeout or error in remote search", "unknown RID or remote search expired",
              "unknown error in remote submission check","timeout in remote submission check function",
              "remote results downloaded and processed", "remote download or post-processing failed", 
              "local search started", "local search finished and processed sucessfully", "local search failed"))

tryCatch({
  
  # ---- Local BLAST ----
  #**********************
  
  #Feedback and Logs
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Start local BLASTn ...\n")}
  newLogs = data.frame(timeStamp = as.integer(Sys.time()), actionId = 1, actionName = "Start local BLASTn")
  
  #Get all files to submit
  if(forceRedo){
    toSubmit = submTable
  } else {
    toSubmit = submTable %>% filter(statusCode == 0)
  }
  
  if(nrow(toSubmit) > 0){
    
    allIds = unique(toSubmit$pipelineId)
    
    for(myId in allIds){
      
      nextSubmit = toSubmit %>% filter(pipelineId == myId)
      
      if(verbose > 0){cat(
        format(Sys.time(), "%H:%M:%S  "), 
        sprintf("- (%i/%i) Process pipelineId %i \n", which(myId == allIds), 
                length(allIds), nextSubmit$pipelineId[1]))}
      
      #Check if the files to blast exist
      if(all(file.exists(paste0(nextSubmit$folder, nextSubmit$fastaFile)))){
        
        #Run blastn first round on on all files to blast
        for(i in 1:nrow(nextSubmit)){
          
          #Update the blastSubmissions table
          myConn = dbConnect(SQLite(), database, synchronous = NULL)
          sqliteSetBusyHandler(myConn, 30000)
          q = dbSendStatement(
            myConn, 
            "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode  = ?, statusMessage = ? WHERE submId = ?",
            params = list(runId, as.integer(Sys.time()), 12, statusCodes$message[12], nextSubmit$submId[i]))
          dbClearResult(q)
          dbDisconnect(myConn)
          
          #Feedback and Logs
          newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                        sprintf("Start BLASTn for pipelineId %i, submId %i", 
                                                nextSubmit$pipelineId[i], nextSubmit$submId[i])))
          if(verbose > 0){cat(paste(rep(" ", 12), collapse = ""),
                              sprintf("BLASTn for submId %i (%i/%i) ... ", 
                                      nextSubmit$submId[i], i, nrow(nextSubmit)))}
          
          # First blast run
          system(sprintf('%s -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -taxidlist %s -num_threads %i -outfmt "%s" | %s > "%s"',
                         blastn, blastArgs$db, paste0(nextSubmit$folder[i], nextSubmit$fastaFile[i]),
                         blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
                         blastArgs$max_hsps, blastArgs$taxidlist, maxCPU, blastArgs$outfmt, zipMethod,
                         paste0(nextSubmit$folder[i], str_replace(nextSubmit$fastaFile[i], ".fasta", ".csv.gz"))))
          
          if(verbose > 0){cat("done\n")}
          
        }
        
        # Second blast run (in case all results identical bit_score)
        blastOut = lapply(
          list.files(nextSubmit$folder[i], full.names = T,
                     pattern = "blastSegmentsClustered\\d+.csv.gz"),
          blast_readOutput, outfmt = blastArgs$outfmt) %>% bind_rows() %>% 
          mutate(geneId = str_extract(qseqid, "^([^_]+)"))
        
        
        
        if(verbose > 0){cat(paste(rep(" ", 12), collapse = ""), "expand BLASTn search ... ")}
        
        #Get segments that have identical blast results (thus might miss some because of 250 lim)
        rerun = blastOut %>% 
          select(segmentId = qseqid, bitscore, geneId, nident, 
                 qlen, length) %>% 
          group_by(segmentId, geneId) %>% 
          filter(bitscore == max(bitscore)) %>% 
          summarise(x = n_distinct(bitscore) == 1, .groups = "drop",
                    identity = nident[1] / qlen[1],
                    cover = length[1] / qlen[1]) %>% 
          filter(x)
        
        #Set these general blast args
        blastArgs$max_target_seqs = 5000 #Set max of 5000 results per target
        blastArgs$qcov_hsp_perc = 100 #Only consider perfect matches
        blastArgs$perc_identity = 100 #Only consider perfect matches
        # blastArgs$taxidlist = sprintf("%s/bact.txids", nextSubmit$folder[1]) #limit to taxa found in first part
        
        
        # write(unique(blastOut$staxids), sprintf("%s/bact.txids", nextSubmit$folder[1]), sep = "\n")
        
        #Blast all perfect matches together
        toFasta = rerun %>% filter(cover == 1, identity == 1)
        
        if(nrow(toFasta) > 0){
          
          myFasta = map_df(1:length(list.files(nextSubmit$folder[i], pattern = "blastSegmentsClustered\\d.csv.gz")), function(x){
            fasta_read(sprintf("%s/blastSegmentsClustered%i.fasta", nextSubmit$folder[1], x),
                       type = "n") %>% filter(id %in% toFasta$segmentId)})
          
          fasta_write(myFasta$seq, sprintf("%s/expand_1.fasta", nextSubmit$folder[1]), 
                      myFasta$id, type = "n")
          
          #Run local blastn on the new database
          system(sprintf('%s -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -num_threads %i -qcov_hsp_perc %.2f -perc_identity %.2f -taxidlist %s -outfmt "%s" | %s > "%s"',
                         blastn, blastArgs$db, sprintf("%s/expand_1.fasta", nextSubmit$folder[1]),
                         blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
                         blastArgs$max_hsps, maxCPU, blastArgs$qcov_hsp_perc,
                         blastArgs$perc_identity, blastArgs$taxidlist, blastArgs$outfmt, zipMethod,
                         sprintf("%s/expand_1.csv.gz", nextSubmit$folder[1])))
        }
        
        #Imperfect ones get blased separately
        toFasta = rerun %>% filter(cover < 1 | identity < 1)
        if(nrow(toFasta) > 0){
          
          myFasta = map_df(1:length(list.files(nextSubmit$folder[i], pattern = "blastSegmentsClustered\\d+.csv.gz")), function(x){
            fasta_read(sprintf("%s/blastSegmentsClustered%i.fasta", nextSubmit$folder[1], x),
                       type = "n") %>% filter(id %in% toFasta$segmentId)})
          
          fasta_write(myFasta$seq, sprintf("%s/expand_2.fasta", nextSubmit$folder[1]), 
                      myFasta$id, type = "n")
          
          #Set the settings to the min cover and identity of imperfect ones
          blastArgs$qcov_hsp_perc = floor(min(toFasta$cover) * 10000) / 10000
          blastArgs$perc_identity = floor(min(toFasta$identity) * 10000) / 10000
          
          #Run local blastn on the new database
          system(sprintf('%s -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -num_threads %i -qcov_hsp_perc %.2f -perc_identity %.2f -taxidlist %s -outfmt "%s" | %s > "%s"',
                         blastn, blastArgs$db, sprintf("%s/expand_2.fasta", nextSubmit$folder[1]),
                         blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
                         blastArgs$max_hsps, maxCPU, blastArgs$qcov_hsp_perc,
                         blastArgs$perc_identity, blastArgs$taxidlist, blastArgs$outfmt, zipMethod,
                         sprintf("%s/expand_2.csv.gz", nextSubmit$folder[1])))
        }
        
        if(verbose > 0){cat("done\n")}
        
        #Update MGS2AMR DB
        #------------------
        
        #Update the DB
        myConn = dbConnect(SQLite(), database, synchronous = NULL)
        sqliteSetBusyHandler(myConn, 30000)
        q = dbSendStatement(
          myConn, 
          sprintf("UPDATE blastSubmissions SET runId = %i, timeStamp = %i, statusCode = %i, statusMessage = \"%s\" WHERE submId IN (%s)",
                  runId, as.integer(Sys.time()), 13, statusCodes$message[13], paste(nextSubmit$submId, collapse = ",")))
        dbClearResult(q)
        
        q = dbSendStatement(
          myConn, 
          "UPDATE pipeline SET statusCode = 4, statusMessage = 'Finished blast', modifiedTimestamp = ?
        WHERE pipelineId = ?",
          params = list(as.character(Sys.time()), nextSubmit$pipelineId[1]))
        dbClearResult(q)
        
        dbDisconnect(myConn)
        
        #Touch the pipelineId to show that the action was completed
        system(paste0("touch ", nextSubmit$folder[1], "pipelineId"))
        
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                      sprintf("Finished BLASTn for submId %s",
                                              paste(nextSubmit$submId, collapse = ","))))
        
      } else {
        
        if(verbose > 0){cat("failed. Files missing\n")}
        
        #Update MGS2AMR DB
        #------------------
        
        #Update the DB
        myConn = dbConnect(SQLite(), database, synchronous = NULL)
        sqliteSetBusyHandler(myConn, 30000)
        q = dbSendStatement(
          myConn, 
          sprintf("UPDATE blastSubmissions SET runId = %i, timeStamp = %i, statusCode = %i, statusMessage = \"%s\" WHERE submId IN (%s)",
                  runId, as.integer(Sys.time()), 14, statusCodes$message[14], paste(nextSubmit$submId, collapse = ",")))
        dbClearResult(q)
        
        dbDisconnect(myConn)
        
        #Touch the pipelineId to show that the action was completed
        system(paste0("touch ", nextSubmit$folder[1], "pipelineId"))
        
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                      sprintf("Failed BLASTn for submId %s",
                                              paste(nextSubmit$submId, collapse = ","))))
        
      }
      
      
    }
  } else {
    
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, "There were no new files to BLAST"))
	  if(verbose > 0){cat(" There were no new files to BLAST\n")}
    
  }
  
  #Feedback and Logs
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Finished local BLASTn"))
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Finished local BLASTn\n")}
  
}, 
finally = {
  
  myConn = dbConnect(SQLite(), database, synchronous = NULL)
  sqliteSetBusyHandler(myConn, 30000)
  
  #Submit the logs, even in case of error so we know where to resume
  newLogs$runId = runId
  newLogs$tool = "localBlast.R"
  newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
  
  q = dbSendStatement(
    myConn, 
    "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
    params = unname(as.list(newLogs)))
  dbClearResult(q)
  dbDisconnect(myConn)
  
})
