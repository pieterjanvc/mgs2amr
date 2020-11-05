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

statusCodes = data.frame(
  status = c(0, 10:13, 20:25, 30:32, 40:42),
  message = c("awaiting submission","successful remote submission", "failed remote submission - server error", 
              "failed remote submission - POST error", "failed remote submission - input error", "searching remotely", 
              "remote search failed","unknown RID or remote search expired", "remote search completed", 
              "unknown error in remote search", "timeout remote search", "remote results download successfull", 
              "remote results download failed", "remote results already downloaded",
              "local search started", "local search finished and processed sucessfully", "local search failed"))


tryCatch({
  
  # ---- Remote BLAST - Submissions ----
  #*************************************
  
  #Feedback and Logs
  newLogs = data.frame(timeStamp = as.integer(Sys.time()), actionId = 1, actionName = "Start remote BLASTn")
  if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Start remote BLASTn ...\n")}
  
  
  #Get all files that need to be submitted
  myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
  toSubmit = dbGetQuery(myConn, paste("SELECT * FROM blastSubmissions WHERE statusCode = 0",
                                      prevRunId)) %>% arrange(timeStamp)
  dbDisconnect(myConn) 
  
  #Submit any new jobs
  if(nrow(toSubmit) > 0){
    
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                                  sprintf("Start submitting %i pending BLASTn jobs", nrow(toSubmit))))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                        sprintf("Start submitting %i pending BLASTn jobs ...\n", nrow(toSubmit)))}

    for(i in 1:nrow(toSubmit)){
      
      #Feedback and Logs
      if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "),
                          sprintf(" - %i/%i Submitting submId %s ... ", i, nrow(toSubmit), toSubmit$submId[i]))}
      Sys.sleep(2) #Prevent too fast consecutive submissions (api will block those)

      #Submit to the blast API 
      mySubmissions = blast_submit(paste0(formatPath(toSubmit$folder[i], endWithSlash = T), toSubmit$fastaFile[i]),
                         program = "blastn", megablast = T, database = "nt", 
                         expect = blastArgs$evalue, word_size = blastArgs$word_size, max_num_seq = blastArgs$max_num_seq,
                         entrez_query = entrezQ, url = blastn, verbose = verbose)

      #Update the blastSubmissions table
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, RID = ?, timeStamp = ?, statusCode  = ?, statusMessage = ?
      WHERE submId = ?",
        params = list(toSubmit$runId[i], RID, as.integer(Sys.time()), mySubmissions$code, 
                      mySubmissions$message, toSubmit$submId[i]))
      dbClearResult(q)
      dbDisconnect(myConn) 

      #Feedback and Logs
      if(mySubmissions$statusCode == 1){
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                                      sprintf("Successful remote submission of %s (submId %s)", toSubmit$fastaFile[i], toSubmit$submId[i])))
        if(verbose > 0){cat("done\n")}
      } else {
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                                      sprintf("Failed remote submission of %s (submId %s)", toSubmit$fastaFile[i], toSubmit$submId[i])))
        if(verbose > 0){cat("failed\n")}
      }
      
    }
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, "Finished submitting pending BLASTn jobs"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "Finished submitting all jobs\n")}
    
  } else {
    #Feedback and Logs
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 6, "No new pending BLASTn jobs, skipping"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "No new pending BLASTn jobs, skipping\n")}
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
                          sprintf("Checking the status of %i active searches ... ", nrow(submTable)))}
      
      #Update submission states
      submTable = submTable %>% rowwise() %>% 
        mutate(blast_checkSubmission(RID, verbose = 0) %>% select(statusCode, statusMessage))
      
      #Update the blastSubmissions table in the DB
      myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
      q = dbSendStatement(
        myConn, 
        "UPDATE blastSubmissions SET runId = ?, timeStamp = ?, statusCode  = ?, statusMessage = ?
        WHERE submId = ?",
        params = list(submTable$runId, rep(as.integer(Sys.time()), nrow(submTable)), submTable$statusCode, 
                      submTable$statusMessage, submTable$submId))
      dbClearResult(q)
      dbDisconnect(myConn)
      
      #Feedback and Logs
      newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, "Finished checking active searches"))
      if(verbose > 0){cat("done\n")}
      
      #Get the submissions that are finished and ready for download 
      submTable = submTable %>% filter(statusCode == 23)
      
      if(nrow(submTable) > 0){
        
        #Feedback and Logs
        newLogs = rbind(newLogs, list(as.integer(Sys.time()), 9, "Start downloading completed searches"))
        if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), 
                            sprintf("Start downloading %i completed searches ...\n", nrow(submTable)))}
        
        for(i in 1:nrow(submTable)){
          
          #Feedback and Logs
          if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "),
                              sprintf(" - %i/%i Downloading results for submId %s ... ", 
                                      i, nrow(submTable), submTable$submId[i]))}
          
          #Get results (success = 9, fail = 10)
          getResults = blast_getResults(submTable$RID[i], submTable$folder[i], unzip = F, verbose = 1)
          
          #Update the blastSubmissions table in the DB
          myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
          q = dbSendStatement(myConn, "UPDATE blastSubmissions 
                    SET statusCode = ?, statusMessage = ? WHERE submId = ?", 
                              params = list(getResults$statusCode, getResults$statusMessage, submTable$submId[i]))
          dbClearResult(q)
          dbDisconnect(myConn)
          
          #Feedback and Logs
          if(getResults$statusCode == 30){
            newLogs = rbind(newLogs, list(as.integer(Sys.time()), 10, 
                                          sprintf("Download/processing successful for BLAST results of %s (submId %s)", 
                                                  submTable$fastaFile[i], submTable$submId[i])))
            if(verbose > 0){cat("done\n")}
          } else {
            newLogs = rbind(newLogs, list(as.integer(Sys.time()), 11, 
                                          sprintf("Download/processing failed for BLAST results of %s (submId %s)", 
                                                  submTable$fastaFile[i], submTable$submId[i])))
            if(verbose > 0){cat("failed\n")}
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
    newLogs = rbind(newLogs, list(as.integer(Sys.time()), 16, "No active searches to check, skipping"))
    if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S "), "No active searches to check, skipping\n")}
  }
  
  #Feedback and Logs
  newLogs = rbind(newLogs, list(as.integer(Sys.time()), 17, "Finished remote BLASTn"))
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


