#!/usr/bin/env Rscript

#******************************
# ---- Local BLASTn search ----
#******************************
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))

#Make sure the blast bin folder is in the R env 
#Sys.getenv("PATH")
#If not, add /opt/ncbi-blast-2.10.1+/bin to PATH in 
#Renviron file at /opt/R/4.0.2/lib/R/etc/Renviron


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

maxCPU = parallel::detectCores()

#Fix and integrate later!!
#The Renv needs to know where the BLASTDB is and this is not taken
#using and export by the user for now ...
Sys.setenv(BLASTDB = blastDB)

#Set these general blast args
blastArgs = list(
  db = "nt",
  evalue = "1e-20",
  word_size = 64,
  max_target_seqs = 250,
  max_hsps = 1,
  taxidlist = sprintf("%sdataAndScripts/%s", baseFolder, "bact.txids"),
  outfmt = "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart send bitscore score length pident nident qcovs qcovhsp"
)

#Limit the search for specific pipelineIds if set, else do all
prevRunId = ifelse(length(pipelineId) != 0,
	sprintf("AND runId in (SELECT runId FROM scriptUse WHERE pipelineId IN ('%s'))", 
	        paste(pipelineId, collapse = "','")),
	"")

#Get all blast submissions that need further follow-up
myConn = dbConnect(SQLite(), database)
submTable = dbGetQuery(myConn, paste("SELECT * FROM blastSubmissions WHERE statusCode in (0,1,3,4)",
									  prevRunId)) %>% arrange(timeStamp)
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
  toSubmit = submTable %>% filter(statusCode == 0)
  if(nrow(toSubmit) > 0){
    for(i in 1:nrow(toSubmit)){
      
      #Update the blastSubmissions table
      myConn = dbConnect(SQLite(), database)
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode  = ?, statusMessage = ? WHERE submId = ?",
        params = list(runId, as.integer(Sys.time()), 12, statusCodes$message[12], toSubmit$submId[i]))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                    sprintf("Start BLASTn for pipelineId %i, submId %i", 
                                            toSubmit$pipelineId[i], toSubmit$submId[i])))
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S  "), 
                          sprintf("- (%i/%i) pipelineId %i : submId %i ... ", 
                                  i, nrow(toSubmit), toSubmit$pipelineId[i], toSubmit$submId[i]))}
      
      #Run local blastn
      system(sprintf('%s -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -taxidlist %s -num_threads %i -outfmt "%s" | gzip > "%s"',
                     blastn, blastArgs$db, paste0(toSubmit$folder[i], toSubmit$fastaFile[i]),
                     blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
                     blastArgs$max_hsps, blastArgs$taxidlist, maxCPU, blastArgs$outfmt, 
                     paste0(toSubmit$folder[i], str_replace(toSubmit$fastaFile[i], ".fasta", ".csv.gz"))))
      
      #Update the DB
      myConn = dbConnect(SQLite(), database)
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode = ?, statusMessage = ?
      WHERE submId = ?",
        params = list(runId, as.integer(Sys.time()), 13, statusCodes$message[13], toSubmit$submId[i]))
      dbClearResult(q)
      
      q = dbSendStatement(
        myConn, 
        "UPDATE pipeline SET statusCode = 4, statusMessage = 'Finished blast', modifiedTimestamp = ?
        WHERE pipelineId = ?",
        params = list(as.character(Sys.time()), toSubmit$pipelineId[i]))
      dbClearResult(q)
      
      dbDisconnect(myConn)
      
      #Touch the pipelineId to show that the action was completed
      system(paste0("touch ", toSubmit$folder[i], "pipelineId"))
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                    sprintf("Finished BLASTn for submId %s",toSubmit$submId[i])))
      if(verbose > 0){cat("done\n")}
      
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
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
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
