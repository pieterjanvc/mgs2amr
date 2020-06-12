library(dplyr)
library(stringr)
library(tidyr)

message("WARNING REPLACE mnt")

#Variables from script
baseFolder = "/mnt/d/Documents/wslData/meta2amr"
inputFile = "/mnt/d/Documents/wslData/test/input.csv"
outputFile = "/mnt/d/Documents/wslData/test/mixedSample1.fastq.gz"
bbmap = "/opt/bbmap"

# baseFolder = "D:/Documents/wslData/meta2amr"
# inputFile = "D:/Documents/wslData/test/input.csv"
# outputFile = "D:/Documents/wslData/test/mixedSample1.fastq.gz"

#Name of the temp folder
tempName = paste0("mixedSample_", as.integer(Sys.time()))

if(file.exists(paste0(baseFolder, "/dataAndScripts/readCounts.csv"))){
  readCounts = read.csv(paste0(baseFolder, "/dataAndScripts/readCounts.csv"),
                        colClasses = c("character", "character", "integer"))
} else {
  readCounts = data.frame(fileName = character(), fileSize = integer(), modDate = character(), readCount = integer())
}

#Read the input file and remove unwanted whitespace
files = read.csv(inputFile, header = T, 
                 colClasses = c("character", "character", "numeric", "character", "character")) %>% 
  mutate(across(where(is.character), function(x) str_trim(x)))
  # mutate(across(where(is.character), function(x) str_trim(x) %>% str_replace("/mnt/d/", "D:/")))


# ---- Check the input file ----
#*******************************

#Check that sum of RA = 1
sumRA = ifelse(sum(files$relativeAbundance) != 1, "*** The sum of relative abundances is not 1", "")

#Type combo check
totalI = sum(str_detect(files$type, "i|I"))
totalB = sum(str_detect(files$type, "b|B"))
isoVsBack = ""
if((totalI < 2 & totalB == 0) | (totalI == 0 & totalB > 0 | (totalB > 1))){
  isoVsBack = "*** Incorrect combination of samples. Choose any of the following:
- Two or more isolate(I) files\n- One background(B) and one or more isolate(I) files"
}


files = files %>% 
  mutate(id = 1:n(), 
         sampleName = ifelse(sampleName == "", paste0("sample", 1:n()), sampleName)) %>% 
  pivot_longer(c(readFile, readFile2), values_to = "filePath") %>% select(id, filePath) %>% 
  filter(filePath != "") %>% 
  mutate(modDate = file.info(filePath)$mtime %>% as.character(), 
         correctType = str_detect(filePath, "\\.fastq\\.gz|\\.fastq")
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


# ---- Get the read counts ----
#*******************************

#Add readcounts from previous runs
files = files %>% mutate(fileSize = file.info(filePath)$size, 
                         fileName = str_extract(filePath, "[^/]+$"))
files = files %>% left_join(readCounts, by = c("fileName", "modDate", "fileSize"))


#If no read counts yet, count them
newFiles = files %>% filter(is.na(readCount)) %>% pull(id) %>% unique()
for(myId in newFiles){
  myFile = files %>% filter(id == myId)
  
  print(paste("Count", myFile$filePath[1]))
  #Count the lines in the file (4 lines = 1 read)
  nReads = system(sprintf("zcat %s | wc -l", myFile$filePath[1]), intern = T) %>%
    as.integer()
  #If painr-end file, total reads id double from counted in one
  nReads = ifelse(nrow(myFile) == 2, nReads / 2, nReads / 4)
  
  # #In case there are two files, interleave them and save in temp
  # if(nrow(myFile) == 2){
  #   reformatOut = system(sprintf("reformat.sh in1=%s in2=%s out=%s/temp/%s/tempFile%s.fastq.gz 2>&1", 
  #                  myFile$filePath[1], myFile$filePath[1], baseFolder, sampleName, myId), 
  #          intern = T)
  #   reformatOut = str_match(reformatOut, "Output:\\s+(\\d+)\\s*reads")[,2]
  #   nreads = reformatOut[!is.na(reformatOut)] %>% as.integer()
  # } else {
  #   nreads = system(sprintf("zcat %s | wc -l", myFile$filePath[1]), intern = T) %>% 
  #     as.integer() / 2
  # }
  
  files[files$id == myId, "readCount"] = nReads
}

#Update the readCounts file
readCounts = rbind(readCounts, files %>% filter(id %in% newFiles) %>% 
                     select(fileName, modDate, fileSize, readCount))
write.csv(readCounts, paste0(baseFolder, "/dataAndScripts/readCounts.csv"), row.names = F)


# ---- Calculate the reads needed for the correct RA ----
#********************************************************

# newCounts = files %>% filter(is.na(readCount)) %>% group_by(id) %>%
#   summarise(filePath = filePath[1], nFiles = n())
# 
# newCounts$readCount = sapply(newCounts$filePath, function(filePath){
# 
#   # rCount = as.integer(system(sprintf("zcat %s | wc -l", filePath), intern = T))
#   rCount = 1
# 
#   if(!is.na(rCount)){
#     return(rCount)
#   } else {
#     stop(paste("The number of reads in", filePath, "can not be counted. Are you sure it is the correct format?"))
#   }
# })
# 
# file.info("D:/Documents/wslData/test/trimmed_read1_SRR4025842.fastq.gz")$size
# 
# files = files %>% left_join(newCounts %>% select(-filePath), by = "id") %>%
#   mutate(readCount = ifelse(is.na(readCount.x), readCount.y, readCount.x)) %>%
#   select(-readCount.x, -readCount.y)
# 
# readCounts = rbind(readCounts, files %>% filter(id %in% newCounts$id) %>% select(filePath, modDate, readCount))
# 
# files = rbind(files %>% filter(!is.na(readCount)), newCounts)
# readCounts = rbind(readCounts, newCounts)
# 
# write.csv(readCounts, paste0(baseFolder, "/dataAndScripts/readCounts.csv"), row.names = F)


# #----------------------------------------
# shell("x=6")
# 
# 
# reformatOut = system(sprintf("/opt/bbmap/reformat.sh in1=%s in2=%s out=%s/test/tempFile%s.fastq.gz 2>&1", 
#                              "/mnt/d/Documents/temp/test/trimmed_read1_SRR4025842.fastq.gz", 
#                              "/mnt/d/Documents/temp/test/trimmed_read2_SRR4025842.fastq.gz", 
#                              "/mnt/d/Documents/temp", 1), 
#                      intern = T)
# 
# 
# reformatOut = system(sprintf("reformat.sh in1=%s in2=%s out=%s/temp/%s/tempFile%s.fastq.gz 2>&1", 
#                              "/data/aplab/ARG_PJ/data/test/trimmed_read1_SRR4025842.fastq.gz", 
#                              "/data/aplab/ARG_PJ/data/test/trimmed_read2_SRR4025842.fastq.gz", 
#                              "/data/aplab/ARG_PJ/data/test", "testSample", 1), 
#                      intern = T)
# 
# reformatOut = str_match(reformatOut, "Output:\\s+(\\d+)\\s*reads")[,2]
# reformatOut[!is.na(reformatOut)] %>% as.integer()
# 663986
# 
# system("/opt/bbmap/reformat.sh", intern = T)
# 
# system(sprintf("zcat %s | wc -l", "/data/aplab/ARG_PJ/data/test/trimmed_read1_SRR4025842.fastq.gz"), intern = T) 