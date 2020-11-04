#!/usr/bin/env Rscript

#******************************
# ---- Local BLASTn search ----
#******************************
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))


# ---- Inputs ----
#*****************

args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file
baseFolder = formatPath(args[[1]], endWithSlash = T)
runId = as.integer(args[[2]])
verbose = args[[3]]
blastn = args[[4]]
blastDB = args[[5]]
prevRunId = as.numeric(strsplit(args[[6]], " "))

#Set these general blast args
blastArgs = list(
  db = "nt",
  evalue = "1e-10",
  word_size = 64,
  max_target_seqs = 50,
  taxidlist = sprintf("%sdataAndScripts/%s", baseFolder, "bact.txids")
)

#Limit the searched for specific runIds if set, else do all
prevRunId = ifelse(!is.na(prevRunId),
	sprintf("AND runId in ('%s')", paste(prevRunId, collapse = "','")),
	"")

#Get all blast submissions that need further follow-up
myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
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

if(system(paste0("if [ -z `command -v ", blastn, "` ]; then echo T; else echo F; fi"), intern = T)){
  stop("The local blastn module cannot be located.\n ",
       "Please check the settings file and update the localBlastBlastn to the blastn module if needed.\n ",
       "Or run remoteBlast.sh if you don't have blast installed (or the resources)")
}

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
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode  = ?, statusMessage = ? WHERE submId = ?",
        params = list(runId, as.integer(Sys.time()), 12, statusCodes$message[12], toSubmit$submId[i]))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                    sprintf("Start BLASTn for submId %i", toSubmit$submId[i])))
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S  "), 
                          sprintf("- Progress %i/%i Blastn for submId %i ... ", 
                                  i, toSubmit$submId[i], nrow(toSubmit)))}
      
      #Run local blastn
      system(sprintf('%s -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -taxidlist %s -outfmt 15 | gzip > "%s"',
                     blastn, blastDB, paste0(toSubmit$folder[i], toSubmit$fastaFile[i]),
                     blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, blastArgs$taxidlist,
                     paste0(toSubmit$folder[i], str_replace(toSubmit$fastaFile[i], ".fasta", ".json.gz"))))
      
      #Update the blastSubmissions table
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode = ?, statusMessage = ?
      WHERE submId = ?",
        params = list(runId, as.integer(Sys.time()), 13, statusCodes$message[13], toSubmit$submId[i]))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      #Add the runId to the temp folder
      write(runId, paste0(toSubmit$folder[i], "runId"))
      
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
