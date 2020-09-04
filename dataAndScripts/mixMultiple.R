#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Load packages
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))

options(digits = 10)
sTime = Sys.time()

# ---- Setup variables and directories ----
#******************************************

#Function to ensure that paths always/never end with '/'
formatPath = function(path, endWithSlash = F){
  
  if(str_detect(path, "\\/$") & endWithSlash){
    return(path)
  } else if(endWithSlash == F){
    return(str_remove(path, "\\/$"))
  } else {
    return(paste0(path, "/"))
  }
  
}

#Variables from mixMultiple.sh script
baseFolder = formatPath(as.character(args[[1]]))
inputFile = as.character(args[[2]])
outputFile = as.character(args[[3]])
readLimit = as.integer(args[[4]])
if(is.na(readLimit)){
  stop("The read limit specified in -m argument needs to be a valid, positive integer")
}
metaData = as.logical(args[[5]])
verbose = as.logical(args[[6]])
tempFolder = formatPath(as.character(args[[7]]))

#Grab the location of the reformat script from the settings file
reformatScript = system(sprintf("grep -oP \"reformat\\s*=\\s*\\K([^\\s]+)\" %s/settings.txt", 
                                baseFolder), intern = T)

#Create temp folder
tempName = paste0("mixedSample_", as.integer(Sys.time()))
tempFolder = paste0(tempFolder, "/", tempName)
dir.create(tempFolder)

if(file.exists(paste0(baseFolder, "/dataAndScripts/readCounts.csv"))){
  readCounts = read.csv(paste0(baseFolder, "/dataAndScripts/readCounts.csv"),
                        colClasses = c("character", "character", "integer", "integer"))
} else {
  readCounts = data.frame(fileName = character(), modDate = character(), 
                          fileSize = integer(), readCount = integer())
}


# ---- Check the input file ----
#*******************************

if(verbose){
  cat(" Check the input file for errors ...")
}

#Read the input file and remove unwanted whitespace
files = read_csv(inputFile, col_names = T, col_types = cols()) %>%  
  mutate(across(where(is.character), function(x) str_trim(x)))

#Check that sum of RA = 1
sumRA = ifelse(sum(files$relativeAbundance) != 1, "*** The sum of relative abundances is not 1", "")

#Type combo check
totalI = sum(str_detect(files$type, "i|I"))
totalB = sum(str_detect(files$type, "b|B"))
isoVsBack = ""
if((totalI < 2 & totalB == 0) | 
   (totalI == 0 & totalB > 0 | 
    (totalB > 1) | 
    sum(totalI, totalB) != nrow(files))){
  isoVsBack = "*** Incorrect combination of samples. Choose any of the following:
- Two or more isolate(I) files\n- One background(B) and one or more isolate(I) files"
}


files = files %>% 
  mutate(id = 1:n(), 
         sampleName = ifelse(sampleName == "", paste0("sample", 1:n()), sampleName)) %>% 
  pivot_longer(c(readFile, readFile2), values_to = "filePath") %>% 
  select(id, type, relativeAbundance, filePath, sampleName) %>% 
  filter(filePath != "") %>% 
  mutate(modDate = file.info(filePath)$mtime %>% as.character(), 
         correctType = str_detect(filePath, "\\.fastq\\.gz$|\\.fastq$")
         )



#Check if files are unique
uniqueFiles = unique(files$filePath)
allFiles = files$filePath

uniqueFiles = ifelse(length(uniqueFiles) != length(allFiles), 
                     paste0("*** The following file names are duplicated:\n", 
                            paste(names(table(allFiles))[table(allFiles) > 1], collapse = "\n")), 
                     "")

#Get all the missing files
missing = c(files %>% filter(is.na(modDate)) %>% pull(filePath) %>% unique())

missing = ifelse(length(missing) > 0, 
                 paste0("*** The following files are missing\n", paste(missing, collapse = "\n")), 
                 "")

#Check for incorrect file types
incorrectType = files %>% filter(!correctType) %>% pull(filePath) %>% unique()
incorrectType = ifelse(length(incorrectType) > 0,paste0("*** The following files are not in fastq or fastq.gz format\n", 
                       paste(incorrectType, collapse = "\n")), 
                       "")
files = files %>% select(-correctType)

#Paste everything together
errorMessage = c(sumRA, isoVsBack, uniqueFiles, missing, incorrectType)
errorMessage = paste(errorMessage[errorMessage != ""], collapse = "\n\n")


if(errorMessage != ""){
  stop(paste0("Incorrect input file\n\n", errorMessage))
}

if(verbose){
  cat(" none found\n")
}

# ---- Get the read counts ----
#*******************************

#Add readcounts from previous runs
files = files %>% mutate(fileSize = file.info(filePath)$size, 
                         fileName = str_extract(filePath, "[^/]+$"))
files = files %>% left_join(readCounts, by = c("fileName", "modDate", "fileSize"))

knownFiles = files %>% filter(!is.na(readCount)) %>% pull(fileName) %>% unique()

if(verbose & length(knownFiles) > 0){
  cat(" Use previous read-count for", knownFiles, "\n")
}

#If no read counts yet, count them
newFiles = files %>% filter(is.na(readCount)) %>% pull(id) %>% unique()

if(length(newFiles) > 0){
  for(myId in newFiles){
    
    myFile = files %>% filter(id == myId)
    
    if(verbose){
      cat(" Counting number of reads in",myFile$fileName, "... ")
    }

    #Count the lines in the file (4 lines = 1 read)
    nReads = system(sprintf("zcat %s | wc -l", myFile$filePath[1]), intern = T) %>%
      as.integer()
    #If pair-end file, total reads id double from counted in one
    nReads = ifelse(nrow(myFile) == 2, nReads / 2, nReads / 4)
    
    
    files[files$id == myId, "readCount"] = nReads
    
    if(verbose){
      cat(nReads, "\n")
    }
  }
  
  #Update the readCounts file
  readCounts = rbind(readCounts, files %>% filter(id %in% newFiles) %>% 
                       select(fileName, modDate, fileSize, readCount))
  write.csv(readCounts, paste0(baseFolder, "/dataAndScripts/readCounts.csv"), row.names = F)
}


# ---- Calculate the reads needed for the correct RA ----
#********************************************************

if(verbose){
  cat(" Calculate the number of reads needed from each file ... ")
}

raData = files %>% group_by(id, type, relativeAbundance, readCount) %>% 
  summarise(.groups = 'drop')

#Get min reads per % 
rpp = min(raData$readCount / (raData$relativeAbundance *  100))

#Calculate the total number of reads
totalReads = sum(raData$relativeAbundance * 100 * rpp)

#If a total is set, adjust the rpp
nReadsM = raData %>% filter(type == "B" | type == "b") %>% pull(readCount)
readLim = ifelse(readLimit == 0 & length(nReadsM) != 0, nReadsM, readLimit)
if(readLim != 0){
  rpp = rpp * readLim / totalReads
}

#Caluclate the times each input file is needed
raData = raData %>% mutate(readsNeeded = relativeAbundance * 100 * rpp,
                           fileNeeded = readsNeeded / readCount)

if(verbose){
  cat("done\n")
}


# ---- Filter and merge the files ----
#*************************************
toMerge = raData %>% left_join(files %>% select(id, filePath), by = "id") %>% 
  group_by(id, fileNeeded) %>% 
  summarise(file1 = filePath[1], file2 = ifelse(is.na(filePath[2]), "", filePath[2]), .groups = 'drop')

for(i in 1:nrow(toMerge)){
  
  if(verbose){
    cat(" Extracting reads from", files %>% filter(id == toMerge$id[i]) %>% pull(fileName), "... ")
  }
  
  fullFile = floor(toMerge$fileNeeded[i])
  partialFile = toMerge$fileNeeded[i] - floor(toMerge$fileNeeded[i])
  
  #Generate a full copy of the file with different ID if file needed more than once
  if(fullFile > 0 & toMerge$fileNeeded[i] != 1.0){
    for(j in 1:fullFile){
      system(sprintf(
        "%s in1=%s in2=%s out=stdout.fastq 2>/dev/null | awk 'NR %% 4 == 1{sub(/@/,\"@%i_\",$0);print;next}\
        NR %% 2 == 1{print \"+\";next}{print}' | gzip -c > %s/tempFile%i_full%i.fastq.gz",
        reformatScript, toMerge$file1[i], toMerge$file2[i],
        i, tempFolder, i, j), intern = F)
    }
  } else if(toMerge$fileNeeded[i] == 1.0){
    system(sprintf(
      "%s in1=%s in2=%s out=%s/tempFile%i_partial.fastq.gz 2>/dev/null",
      reformatScript, toMerge$file1[i], toMerge$file2[i],
      tempFolder, i), intern = T)
  }
  
  #Filter the fraction of reads needed 
  if(partialFile != 0){
    system(sprintf(
      "%s --samplerate=%0.10f in1=%s in2=%s out=%s/tempFile%i_partial.fastq.gz 2>/dev/null",
      reformatScript, partialFile, toMerge$file1[i], toMerge$file2[i],
      tempFolder, i), intern = T)
  }
  
  if(verbose){
    cat("done\n")
  }
  
}

#Merge the temp files into the final one
if(verbose){
  cat(" Merge all reads together and write final file ... ")
}

system(paste0("cat ", tempFolder, "/*.fastq.gz > ", outputFile))

#Write the meta data as JSON (if requested)
if(metaData){
  raData$readsNeeded = as.integer(raData$readsNeeded)
  
  metaData = list(
    timestamp = Sys.time() %>% as.character(),
    inputFile = inputFile,
    outputFile = outputFile,
    readLimit = readLimit,
    totalReads = sum(raData$readsNeeded),
    fileData = raData %>% left_join(files %>% select(-type, -relativeAbundance, -readCount), by = "id") %>% 
      group_by(across(c(-filePath, -modDate, -fileSize, -fileName))) %>%
      summarise(fileName1 = fileName[1], 
                filePath1 = filePath[1], 
                fileName2 = ifelse(is.na(fileName[2]), "", fileName[2]), 
                filePath2 = ifelse(is.na(filePath[2]), "", filePath[2]), .groups = "drop")
  )

  write_json(metaData, paste0(str_extract(outputFile, ".*(?=\\.fastq\\.gz$)"), "_metaData.json"), pretty = T)
  
}

#Remove temp files
system(paste("rm -r", tempFolder))

if(verbose){
  cat("done\n\n Total time to run script:", 
      round(difftime(Sys.time(), sTime, units = "mins"),2), "minutes",
      "\n")
  
}
