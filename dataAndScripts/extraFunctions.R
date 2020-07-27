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

#Submit to BLAST API
submitBlast = function(fastaPath, program = "blastn", database = "nt", expect = 1e-10, word_size = 64, max_num_seq = 50, 
                       entrez_query = NULL, entrez_query_preset_excl = NULL, eq_menu = NULL, eq_op = NULL, exclude_seq_uncult = NULL,
                       verbose = T){
  
  #Build the URL with settings
  blastURL = parse_url("https://blast.ncbi.nlm.nih.gov/Blast.cgi")
  blastURL$query = list(
    PROGRAM = program, 
    DATABASE = database, 
    EXPECT = expect, 
    WORD_SIZE = word_size, 
    MAX_NUM_SEQ = max_num_seq, 
    ENTREZ_QUERY = entrez_query, 
    ENTREZ_QUERY_PRESET_EXCL = entrez_query_preset_excl, 
    EQ_MENU = eq_menu, 
    EQ_OP = eq_op, 
    EXCLUDE_SEQ_UNCULT = exclude_seq_uncult
  )
  
  submitToNCBI = POST(build_url(blastURL),
                      body = list(
                        QUERY = paste(readLines(fastaPath), collapse = "\n")
                      ), 
                      encode = "multipart")
  
  # Check if successful
  if(submitToNCBI$status_code == 200){
    # Extract the search ID 
    RID = content(submitToNCBI) %>% as.character() %>% str_extract("(?<=RID\\s=\\s)\\w+")
    if(!is.na(RID)){
      return(RID)
    } else {
      message("The submission was returned with no RID. Please check all arguments to make sure they are correct")
      return(NA)
    }
  } else {
    message(paste("Submission failed with status code", submitToNCBI$status_code))
    return(NA)
  }
  
}

#Check the status of a submission
checkBlastSubmission = function(RID, keepChecking = T, checkInterval = 30, timeOut = 3600, verbose = T){
  
  searching = T
  startTime = Sys.time()
  result = F
  
  while(searching){
    mySearch = GET(paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=", RID))
    mySearch = content(mySearch, encoding = "latin1") %>% as.character() %>% str_extract("(?<=Status=)\\w+")
    
    # Evaluate search progress
    if(mySearch == "WAITING"){
      if(verbose){message(sprintf("Searching ... %.1f minutes since start", difftime(Sys.time(), startTime, "mins")))}
    } else if(mySearch == "FAILED"){
      if(verbose){message("Search failed - stop")}
      searching = F
    } else if(mySearch == "UNKNOWN"){
      if(verbose){message("Unknown RID or search expired - stop")}
      searching = F
    } else if(mySearch == "READY"){
      if(verbose){message("Search completed!")}
      result = T
      searching = F
    } else {
      if(verbose){message("Unknown error  - stop")}
      searching = F
    }
    
    if(difftime(Sys.time(), startTime, "sec") >= timeOut){
      if(verbose){message("Timeout reached - stop checking")}
      searching = F
    }
    
    if(!keepChecking | !searching){
      searching = F
    } else {
      Sys.sleep(checkInterval)
    }
    
  }
  
  #Return TRUE if search successfully completed, else FALSE
  return(result)
}

#Get the results for an RID and save as JSON
getBlastResults = function(RID, outputFolder, unzip = T, mergeMultiple = T, omitSeqData = F, 
                           returnJSON = T, verbose = T){
  
  outputFolder = formatPath(outputFolder)
  
  if(checkBlastSubmission(RID, verbose = F, keepChecking = F)){ #Check if RID exists
    
    #Download...
    if(verbose){message("Downloading results ... ", appendLF = F)}
    myResult = GET(sprintf("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=JSON2&RID=%s", RID),
                   write_disk(paste0(outputFolder, RID, ".gz"), overwrite = F))
    if(verbose){message("done")}
    
    if(unzip){ #Unzip if set T
      
      if(verbose){message("Unzip data ... ", appendLF = F)}
      untar(paste0(outputFolder, RID, ".gz"), exdir = paste0(outputFolder, RID))
      if(verbose){message("done")}

      if(mergeMultiple){ #Merge multiple files into 1 JSON
        if(verbose){message("Merging and updating results ... ", appendLF = F)}
        file.remove(paste0(paste0(outputFolder, RID, "/"), RID, ".json"))        
        
        myJSON = sapply(list.files(paste0(outputFolder, RID, "/"), full.names = T), readLines) %>% unlist() %>% unname()
        
        if(omitSeqData){
          #Remove seq and other large attributes (not needed)
          findLines = str_detect(myJSON, "^\\s+(\"qseq\":)|(\"hseq\":)|(\"midline\":)")
          myJSON = myJSON[!findLines]
          #Fix trailing comma on last one before removed
          findLines = str_detect(myJSON, "^\\s+\"gaps")
          myJSON[findLines] = str_remove(myJSON[findLines], ",$")
        }
        
        #Update the top names of the JSON files
        findLines = str_detect(myJSON, "^\\s+\"BlastOutput2")
        query_titles = str_detect(myJSON, "^\\s+\"query_title")
        if(sum(query_titles) > 0){
          query_titles = str_extract(myJSON[query_titles], "(?<=\")[^\"]+(?=\",)")
        } else {
          query_titles = str_detect(myJSON, "^\\s+\"query_id")
          query_titles = str_extract(myJSON[query_titles], "(?<=\")[^\"]+(?=\",)")
        }
        myJSON[findLines] = paste0("  \"read", query_titles, "\": {")
        #Remove unwanted {
        findLines = str_detect(myJSON, "^\\{")
        myJSON = myJSON[!c(F, findLines[-1])]
        #Remove unwanted }
        findLines = str_detect(myJSON, "^\\}")
        myJSON = myJSON[!c(findLines[-length(findLines)], F)]
        #Add , to separate different files (ignore last)
        findLines = str_detect(myJSON, "^  \\}")
        findLines[length(findLines)-1] = F
        myJSON[findLines] = "  },"
        
        if(verbose){message("done")}
        
        if(verbose){message("Write data and clean-up ... ", appendLF = F)}
        #Save
        write(myJSON, paste0(outputFolder, RID, ".json"))
        
        #Remove temp folder
        unlink(paste0(outputFolder, RID, "/*"), recursive = T)
        unlink(paste0(outputFolder, RID), recursive = T)
        
        if(returnJSON){
          if(verbose){message("done")}
          return(read_json(paste0(outputFolder, RID, ".json")))
        } else {
          if(verbose){message("done")}
          return(T)
        }
        
      } else {
        return(T)
      }
      
    } else {
      return(T)
    }
  } else {
    message("No result found. Check the RID")
    return(F)
  }
  
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

#Make cut-offs based on steepest tangent in curve (maybe not best if very different RA)
cutOff = function(numbers){
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique()
}


writeFasta = function(id, seq, filePath, extraInfo, type, checkData = T){
  
  # https://en.wikipedia.org/wiki/FASTA_format
  
  #INPUT CHECKS
  if(checkData){
    checkInput = table(id) > 1
    if(any(checkInput)){
      stop(paste("The IDs must be unique. \n Violating IDs:", paste(names(checkInput[checkInput]), collapse = ", ")))
    }
    
    checkInput = str_detect(id, "\\s")
    if(any(checkInput)){
      stop(paste("The IDs cannot contain spaces.\n Use the extraInfo argument",
                 "if additional text needs to be in the header.\n Violating IDs:", 
                 paste(id[checkInput], collapse = ", ")))
    }
    
    if(tolower(type) == "dna"){
      checkInput = str_detect(seq, "[^ACGTURYKMSWBDHVN-]")
    } else if(tolower(type) %in% c("prot", "protein")) {
      checkInput = str_detect(seq, "[^A-Z*-]")
    } else {
      stop("type needs to be either 'dna' or 'prot'")
    }
    
    if(any(checkInput)){
      stop(paste("The following IDs have seqences that are not valid", paste(id[checkInput], collapse = ", ")))
    }
    
    if(!str_detect(filePath, "\\.(fasta|fa|fna|ffn|faa|frn)$")){
      warning("The file defined had not a common fasta extension")
    }
  }

  #Generate the fasta
  myFasta = data.frame(
    id = paste0(">", id, if(!missing(extraInfo)){paste0(" ", extraInfo)}), 
    seq = seq) %>% select(id, seq) %>% pivot_longer(cols = everything()) %>% select(value)
  
  #Write the file
  write.table(myFasta, file = filePath, sep = "\t", 
              quote = F, row.names = F, col.names = F)
  
}

readFasta = function(filePath, checkData = T, type){
  
  #Read fasta file per line
  myFasta = readLines(filePath)
  
  #Find headers
  headers = which(str_detect(myFasta, "^>"))
  
  #Extract sequences
  seq = sapply(1:(length(headers)-1), function(x){
    paste0(myFasta[(headers[x]+1):(headers[x+1]-1)])
  })
  seq = c(seq, paste0(myFasta[(headers[length(headers)]+1):length(myFasta)]))
  headers = myFasta[headers]
  
  #Create data frame
  myFasta = data.frame(id = str_extract(headers, "(?<=>)[^\\s]+"), seq = seq, 
             extraInfo = str_remove(headers, ">[^\\s]+\\s*"))
  
  #Check data if set
  if(checkData){
    checkInput = table(myFasta$id) > 1
    if(any(checkInput)){
      warning(paste("The IDs are not unique. \n Violating IDs:", paste(names(checkInput[checkInput]), collapse = ", ")))
    }
    
    if(tolower(type) == "dna"){
      checkInput = str_detect(myFasta$seq, "[^ACGTURYKMSWBDHVN-]")
    } else if(tolower(type) %in% c("prot", "protein")) {
      checkInput = str_detect(myFasta$seq, "[^A-Z*-]")
    } else {
      stop("type needs to be either 'dna' or 'prot'")
    }
    
    if(any(checkInput)){
      warning(paste("The following IDs have seqences that are not valid\n",
                    "(Make sure you have set the correct 'type' argument (dna or prot))\n",
                    paste(myFasta$id[checkInput], collapse = ", ")))
    }
  }
  
  myFasta
}
