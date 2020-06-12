#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file

nReadsM = as.numeric(args[[1]])
nReadsI = as.numeric(unlist(strsplit(args[[2]], ",")))
percentMixIn = as.numeric(unlist(strsplit(args[[3]], ",")))
limit = as.numeric(args[[4]])

# nReadsM = 28452558
# nReadsI = c(2799742,983158,2523678,1750458,2690716)
# percentMixIn = c(0.1,0.15,0.3,0.2,0.1)
# procedure = "oversample"
# nReadsI / sum(nReadsI)

# nReadsM = 0
# nReadsI = c(200, 180)
# percentMixIn = c(0.4, 0.6)
# limit = 400

#If there is a background metagenome, add this to the calculations
if(nReadsM !=0){
  percentMixIn = c(percentMixIn, 1 - sum(percentMixIn))
  nReadsI = c(nReadsI, nReadsM)
}

#Get min reads per % 
rpp = min(nReadsI / (percentMixIn *  100))

#Calculate the total number of reads
totalReads = sum(percentMixIn * 100 * rpp)

#If a total is set, adjust the rpp
limit = ifelse(limit == 0 & nReadsM != 0, nReadsM, limit)
if(limit != 0){
  rpp = rpp * limit / totalReads
}

#Caluclate the times each input file is needed
result = (percentMixIn * 100 * rpp) / nReadsI

#Sanity check
#result * nReadsI

#Return percentages per file (last is background)
cat(sapply(result, function(x) sprintf("%.6f", x)))
