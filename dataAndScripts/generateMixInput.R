#!/usr/bin/env Rscript

# ---- Generate mixMultiple.sh input file ----
#*********************************************
library(stringr)
args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file

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

tempFolder = formatPath(args[[1]])
bgFolder = formatPath(args[[2]])
bgName = args[[3]]
iFolder = formatPath(args[[4]])
iNames = str_split(args[[5]], ",") %>% unlist()
relAb = str_split(args[[6]], ",") %>% unlist() %>% as.numeric()
fileName = args[[7]]

# #For testing
# tempFolder = "S:/pipelineTemp/"
# bgFolder = "P:/ARG_PJ/aim2/meta2amr/dataAndScripts/testData/"
# bgName = "background"
# iFolder = "P:/ARG_PJ/aim2/meta2amr/dataAndScripts/testData/"
# iNames = c("isolate1", "isolate2")
# relAb = c(0.1,0.3)
# fileName = "test"

#Get the files names (can be 1 or 2 per sample)
bgFiles = if(bgName == ""){NA} else {list.files(bgFolder, pattern = bgName, full.names = T)}
iFiles = list.files(iFolder, pattern = ".fastq.gz", full.names = T)
iFiles = lapply(iNames, function(x) iFiles[str_detect(iFiles, x)])

#Build the csv file
myFile = data.frame(
  type = "I", sampleName = iNames,
  relativeAbundance = relAb,
  readFile = sapply(iFiles, "[", 1), 
  readFile2 = sapply(iFiles, "[", 2)
)

if(bgName != ""){
  myFile[nrow(myFile) + 1,] = list("B", bgName, 1 - sum(relAb), bgFiles[1], bgFiles[2])
}

#Write the csv file
write.csv(myFile, paste0(tempFolder, "/", fileName, ".csv"), row.names = F)
