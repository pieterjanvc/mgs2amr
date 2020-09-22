#!/usr/bin/env Rscript

#**********************
# ---- Blast prep ----
#*********************
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gfaTools))
suppressPackageStartupMessages(library(RSQLite))


# ---- Inputs ----
#*****************
args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file

baseFolder = formatPath(args[[1]], endWithSlash = T)
verbose = args[[2]]

#Get all blast submissions that need further follow-up
myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
submTable = dbGetQuery(myConn, "SELECT * FROM blastSubmissions WHERE statusCode in (1,3,4)")
dbDisconnect(myConn) 

#Status codes
statusCodes = data.frame(
  status = 1:10,
  message = c("successful submission", "failed submission", "searching", 
              "search completed", "timeout or error in search", "unknown RID or search expired",
              "unknown error in submission check","timeout in submission check function",
              "results downloaded and processed", "download or post-processing failed")
)


#Update submission states
submTable = submTable %>% 
  mutate(statusCode = ifelse(statusCode %in% c(1,3), 
                             blast_checkSubmission(RID, verbose = 0)$statusCode,
                             statusCode))

#Download the results for all newly completed searches
submTable = submTable %>%
  mutate(statusCode = 
           ifelse(statusCode == 4,
                  ifelse(blast_getResults(RID, folder, omitSeqData = T, verbose = 1, returnJSON = F),
                         9, 10),
                  statusCode))

#Update status messages
submTable = submTable %>%
  mutate(statusMessage = statusCodes$message[statusCode])

#Update the database (if any changes)
if(nrow(submTable) > 0){
  myConn = dbConnect(SQLite(), paste0(baseFolder, "dataAndScripts/meta2amr.db"))
  q = dbSendStatement(myConn, "UPDATE blastSubmissions 
                    SET statusCode = ?, statusMessage = ? WHERE RID = ?", 
                      params = list(submTable$statusCode, submTable$statusMessage, submTable$RID ))
  dbClearResult(q)
  dbDisconnect(myConn)
}

