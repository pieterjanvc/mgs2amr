#********************
# ---- FUNCTIONS ----
#********************

library(httr)
library(dplyr)
library(jsonlite)
library(stringr)

# Format a path so it always end with or without slash
formatPath = function(path, endWithSlash = T){
  
  if(str_detect(path, "\\/$") & endWithSlash){
    return(path)
  } else if(endWithSlash == F){
    return(str_remove(path, "\\/$"))
  } else {
    return(paste0(path, "/"))
  }
  
}

# Call the BLAST API (NCBI ONLINE)
BLASTapi = function(fastaPath, resultsFolder, waitUntilDone = T, secBetweenCheck = 30){
  
  # Submit the file through API
  submitToNCBI = POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Put&PROGRAM=blastn&DATABASE=nt&LINK_LOC=blasthome",
                      body = list(
                        EXPECT = 1e-10, 
                        WORD_SIZE = 64,
                        MAX_NUM_SEQ=50,
                        ENTREZ_QUERY= "(txid2+%5BORGN%5D)+NOT(environmental+samples%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+txid32644%5Borgn%5D)",
                        ENTREZ_QUERY_PRESET_EXCL= "environmental+samples%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+txid32644%5Borgn%5D&EQ_MENU=bacteria+(taxid%3A2)",
                        EQ_OP= "AND",
                        EXCLUDE_SEQ_UNCULT= "on",
                        QUERY = paste(readLines(fastaPath), collapse = "\n")
                        ), 
                      encode = "multipart")
  
  # Check if successful
  if(submitToNCBI$status_code == 200){
    # Extract the search ID 
    RID = content(submitToNCBI) %>% as.character() %>% str_extract("(?<=RID\\s=\\s)\\w+")
    
    if(waitUntilDone){
      
      # Check NCBI every so often to see if search is completed
      Sys.sleep(30) # start with 30 sec wait (api priority is low and slow!)
      searching = T
      startTime = Sys.time()
      while(searching){
        Sys.sleep(secBetweenCheck)
        mySearch = GET(paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=", RID))
        mySearch = content(mySearch) %>% as.character() %>% str_extract("(?<=Status=)\\w+")
        
        # Evaluate search progress
        if(mySearch == "WAITING"){
          print(sprintf("Searching ... %.1f minutes since submission", difftime(Sys.time(), startTime, "mins")))
        } else if(mySearch == "FAILED"){
          print("Search failed - stop")
          searching = F
        } else if(mySearch == "UNKNOWN"){
          print("Search expired - stop")
          searching = F
        } else if(mySearch == "READY"){
          print("Search done!")
          searching = F
        } else {
          print("Unknown error  - stop")
          searching = F
        }
        
      }
      
      # If successful, retrieve results and save to disk (zip)
      if(mySearch == "READY"){
        mySearch = GET(sprintf("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=JSON2&RID=%s", RID),
                       write_disk(paste0(formatPath(resultsFolder), RID, ".gz"), overwrite = T))
      } else {
        print("Search failed")
      }
    }
    
    
  } else {
    print("Submission to NCBI failed")
  }
  
  return(RID)
}


#Get the coverage based
getCoverage = function(from, to, readLength){
  ranges = data.frame(from = from, to = to) %>% 
    arrange(from, to)
  totalLength = ranges$to[1] - ranges$from[1] + 1
  maxPos = ranges$to[1]
  if(nrow(ranges) > 1){
    for(i in 2:nrow(ranges)){
      
      if(ranges$from[i] <= ranges$to[i-1] & 
         ranges$to[i] > maxPos) {
        totalLength = totalLength + ranges$to[i] - maxPos
      } else if(ranges$from[i] > ranges$to[i-1] & 
                ranges$to[i] > maxPos){
        totalLength = totalLength + ranges$to[i] - max(maxPos, ranges$from[i])
      }
      
      maxPos = max(maxPos, ranges$to[i])
      
    }
  }
  
  return(totalLength / readLength)
}