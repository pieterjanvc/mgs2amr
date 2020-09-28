#!/usr/bin/env Rscript

#******************************
# ---- Remote BLASTn search ----
#******************************
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))


# ---- Inputs ----
#*****************

#Set these general blast args
blastArgs = list(
  db = "nt",
  evalue = "1e-10",
  word_size = 64,
  max_num_seq = 50
)

args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file
baseFolder = formatPath(args[[1]], endWithSlash = T)
runId = as.integer(args[[2]])
verbose = args[[3]]
blastn = args[[4]]
entrezQ = args[[5]]
prevRunId = as.numeric(strsplit(args[[6]], " "))
timeOut = as.integer(args[[7]])
checkFreq = as.integer(args[[8]])

# baseFolder = "./"
# runId = 1
# verbose = 1
# blastn = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# entrezQ = "txid2 [ORGN]"
# prevRunId = NA
# timeOut = 3600
# checkFreq = 120

start_time = as.integer(Sys.time())
#Check if the URL is working (error when not)
test = blast_checkSubmission("test", url = blastn, verbose = 0)

#Limit the searched for specific runIds if set, else do all
prevRunId = ifelse(!is.na(prevRunId),
                   sprintf("AND runId in ('%s')", paste(prevRunId, collapse = "','")),
                   "")

#Status codes
statusCodes = data.frame(
  status = 0:13,
  message = c("awaiting submission","successful remote submission", "failed remote submission", "searching remotely", 
              "remote search completed", "timeout or error in remote search", "unknown RID or remote search expired",
              "unknown error in remote submission check","timeout in remote submission check function",
              "remote results downloaded and processed", "remote download or post-processing failed", 
              "local search started", "local search finished and processed sucessfully", "local search failed"))


tryCatch({
  
  # ---- Remote BLAST - Submissions ----
  #*************************************
  
  #Feedback and Logs
  newLogs = data.frame(timeStamp = as.integer(Sys.time()), actionId = 1, actionName = "Start remote BLASTn")
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Start remote BLASTn ...\n")}
  
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, "Start submitting pending BLASTn jobs"))
  
  #Get all files that need to be submitted
  myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
  toSubmit = dbGetQuery(myConn, paste("SELECT * FROM blastSubmissions WHERE statusCode = 0",
                                      prevRunId)) %>% arrange(timeStamp)
  dbDisconnect(myConn) 
  
  #Submit any new jobs
  if(nrow(toSubmit) > 0){
    
    #Feedback and Logs
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                        sprintf("Start submitting %i pending BLASTn jobs ...\n", nrow(toSubmit)))}
    
    for(i in 1:nrow(toSubmit)){
      
      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "),
                          sprintf(" - Submitting %s (submId %s) ...", toSubmit$fastaFile, toSubmit$submId))}
      Sys.sleep(2) #Prevent too fast consecutive submissions (api will block those)
      
      #Submit to the blast API 
      RID = blast_submit(paste0(formatPath(toSubmit$folder[i], endWithSlash = T), toSubmit$fastaFile[i]),
                         program = "blastn", megablast = T, database = "nt", 
                         expect = blastArgs$evalue, word_size = blastArgs$word_size, max_num_seq = blastArgs$max_num_seq,
                         entrez_query = entrezQ, url = blastn, verbose = verbose)
      #Status depends on (un)successful submission
      statusCode = ifelse(is.na(RID), 2, 1)
      
      #Update the blastSubmissions table
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, RID = ?, timeStamp = ?, statusCode  = ?, statusMessage = ?
      WHERE submId = ?",
        params = list(toSubmit$runId[i], RID, as.integer(Sys.time()), statusCode, 
                      statusCodes$message[statusCode+1], toSubmit$submId[i]))
      dbClearResult(q)
      dbDisconnect(myConn) 
      
      #Feedback and Logs
      if(statusCode == 1){
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                      sprintf("Successful remote submission of %s (submId %s)", toSubmit$fastaFile, toSubmit$submId)))
        if(verbose > 0){"done\n"}
      } else {
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                                      sprintf("Failed remote submission of %s (submId %s)", toSubmit$fastaFile, toSubmit$submId)))
        if(verbose > 0){"failed\n"}
      }
      
    }
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Finished submitting pending BLASTn jobs"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Finished submitting all jobs\n")}
    
  } else {
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, "No pending BLASTn jobs"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "No pending BLASTn jobs\n")}
  }
  
  
  # ---- Remote BLAST - Retrieval of results ----
  #**********************************************
  
  #Get all blast submissions that need further follow-up
  myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
  submTable = dbGetQuery(myConn, paste("SELECT * FROM blastSubmissions WHERE statusCode in (1,3,4)",
                                       prevRunId)) %>% arrange(timeStamp)
  dbDisconnect(myConn)
  
  if(nrow(submTable) > 0){
    
    #Run this loop as long as there are pending searches and timeout not reached
    while(nrow(submTable) > 0 & (as.integer(Sys.time()) - start_time < timeOut)){
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 7, "Start checking active searches"))
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                          sprintf("Checking the status of %i active searches ...", nrow(submTable)))}
      
      #Update submission states
      submTable = submTable %>% rowwise() %>% 
        mutate(statusCode = blast_checkSubmission(RID, verbose = 0)$statusCode)
      
      #Update the blastSubmissions table in the DB
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode  = ?, statusMessage = ?
        WHERE submId = ?",
        params = list(submTable$runId, rep(as.integer(Sys.time()), nrow(submTable)), submTable$statusCode, 
                      statusCodes$message[submTable$statusCode+1], submTable$submId))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, "Finished checking active searches"))
      if(verbose > 0){"done\n"}
      
      #Get the submissions that are finished and ready for download 
      submTable = submTable %>% filter(statusCode == 4)
      
      if(nrow(submTable) > 0){
        
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Start downloading completed searches"))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                            sprintf("Start downloading %i completed searches ...\n", nrow(submTable)))}
        
        for(i in 1:nrow(submTable)){
          
          #Feedback and Logs
          if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "),
                              sprintf(" - %i/%i Downloading results for submId %s ...", 
                                      i, nrow(submTable), toSubmit$submId))}
          
          #Get results (success = 9, fail = 10)
          statusCode = ifelse(blast_getResults(submTable$RID[i], submTable$folder[i], omitSeqData = T, verbose = 1, returnJSON = F),
                              9, 10)
          
          #Update the blastSubmissions table in the DB
          myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
          q = dbSendStatement(myConn, "UPDATE blastSubmissions 
                    SET statusCode = ?, statusMessage = ? WHERE submId = ?", 
                              params = list(statusCode, statusCodes$message[statusCode+1], submTable$submId[i]))
          dbClearResult(q)
          dbDisconnect(myConn)
          
          #Feedback and Logs
          if(statusCode == 9){
            newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, 
                                          sprintf("Download/processing successful for BLAST results of %s (submId %s)", 
                                                  submTable$fastaFile, submTable$submId)))
            if(verbose > 0){"done\n"}
          } else {
            newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, 
                                          sprintf("Download/processing failed for BLAST results of %s (submId %s)", 
                                                  submTable$fastaFile, submTable$submId)))
            if(verbose > 0){"failed\n"}
          }
        }
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 12, "Finished downloading completed searches"))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Finished downloading completed searches\n")}
        
      } else {
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 13, "No new downloads found, skipping"))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "No new downloads found, skipping\n")}
        
      }
      
      #Get the new status of all submissions
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      submTable = dbGetQuery(myConn, paste("SELECT * FROM blastSubmissions WHERE statusCode in (1,3,4)",
                                           prevRunId)) %>% arrange(timeStamp)
      dbDisconnect(myConn)
      
      #Check if the timeout has not been reached, if so, exit the script
      if(as.integer(Sys.time()) - start_time > timeOut){
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 14, sprintf("Script timeout reached (%i sec)", timeOut)))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), sprintf("Script timeout reached (%i sec)\n", timeOut))}
        
      } else if(nrow(submTable) > 0){
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 15, 
                                      sprintf("Waiting %i sec before checking again for new results", checkFreq)))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                            sprintf("%i results pending: Waiting %i sec before checking again for new results\n", 
                                    nrow(submTable), checkFreq))}
        Sys.sleep(checkFreq)
      }
    }
  } else {
    
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 16, "No active searches to check"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "No active searches to check\n")}
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 17, "Skip downloading as there are no active searches"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Skip downloading as there are no active searches\n")}
  }
  
  #Feedback and Logs
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 18, "Finished remote BLASTn"))
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Finished remote BLASTn\n")}
  
}, 
finally = {
  
  myConn = dbConnect(SQLite(), sprintf("%sdataAndScripts/meta2amr.db", baseFolder))
  
  #Submit the logs, even in case of error so we know where to resume
  newLogs$runId = runId
  newLogs$tool = "remoteBlast.R"
  newLogs = newLogs %>% select(runId,tool,timeStamp,actionId,actionName)
  
  q = dbSendStatement(myConn, "INSERT INTO logs (runId,tool,timeStamp,actionId,actionName) VALUES(?,?,?,?,?)", 
                      params = unname(as.list(newLogs)))
  dbClearResult(q)
  dbDisconnect(myConn)
  
})


