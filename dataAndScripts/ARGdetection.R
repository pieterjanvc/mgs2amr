

library(stringr)
library(dplyr)
library(purrr)
library(igraph)
library(gfaTools)
library(tidyr)
library(jsonlite)

#cat $(find . -name '*.gfa') > masterGFA.gfa

# ---- Merge all MC output ----
#******************************

baseFolder = "D:/Documents/wslData/meta2amr"
tempName = "testOutput_1592520811"
tempName = "EC_KP_MC"

#Read the master GFA file as gfa object
filePath = paste0(baseFolder, "/temp/", tempName, "/masterGFA.gfa")
gfa = readGFA(filePath)

#Read the raw file
myFile = str_split(readLines(filePath), "\t")

#Get all IDs for which the file is not empty
message("REMOVE WSL ON LINUX")
baseFolder = "/mnt/d/Documents/wslData/meta2amr"
geneId = system(paste("wsl find ", paste0(baseFolder, "/temp/", tempName, "/"), 
                    " -name '*.gfa' | xargs wc -l"), intern = T)
geneId = str_match(geneId, "(\\d+).*/(\\d+)/graph\\.gfa$")
geneId = data.frame(geneId = geneId[,3], count = as.integer(geneId[,2])) %>% 
  na.omit() %>% filter(count != 0) %>% pull(geneId)


# #The folder number is the gene ID
# geneId = system(paste0("ls -d ",baseFolder, "/temp/", tempName, "/*/"), intern = T) %>% 
#   str_extract("\\d+(?=/$)") %>% as.integer()
# geneId = geneId[!is.na(geneId)]

#Add the gene ID to the segments
start = which(sapply(myFile, "[[", 1) == "S")
diff = c(start, 1) - c(-2, start)
start = diff > 1
pos = which(diff > 1)

gfa$segments$geneId = rep(geneId, c(pos[-1], length(start)) - pos)

#Add the gene ID to the links
start = which(sapply(myFile, "[[", 1) == "L")
diff = c(start, 1) - c(-2, start)
start = diff > 1
pos = which(diff > 1)

gfa$links$geneId = rep(geneId, c(pos[-1], length(start)) - pos)

#REMOVE LATER
baseFolder = "D:/Documents/wslData/meta2amr"
rm(myFile)

#Save the intermediate gfa file 
saveRDS(gfa, file = paste0(baseFolder, "/temp/", tempName, "/masterGFA.rds"))
gfa = readRDS(paste0(baseFolder, "/temp/", tempName, "/masterGFA.rds"))


# ---- Detect important ARG ----
#*******************************

#Extract the kmercounts
kmerCounts = gfa$segments %>% 
  mutate(start = str_detect(name, "_start$")) %>% 
  rename(readId = name)

#Get the ARG list
argGenes = read.csv(paste0(baseFolder, "/dataAndScripts/argTable.csv"), 
                    colClasses = list(geneId = "character"))

#Detect the most likely genes
genesDetected = kmerCounts %>% filter(start) %>% 
  left_join(argGenes %>% select(geneId, gene, subtype, nBases, name), by = c("geneId" = "geneId")) %>% 
  group_by(geneId, gene, subtype, nBases, name) %>% 
  summarise(length = sum(LN), kmerCount = sum(KC), n = n(), .groups = 'drop') %>% rowwise() %>% 
  mutate(covered = min(1, length / nBases)) %>% 
  group_by(gene) %>% filter(kmerCount == max(kmerCount))

#Make cut-offs based on steepest tangent in curve (maybe not best if very different RA)
cutOff = function(numbers){
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique()
}

genesDetected = genesDetected %>% filter(covered >= cutOff(genesDetected$covered)) 

write.csv(genesDetected, paste0(baseFolder, "/temp/", tempName, "/genesDetected.csv"))


# ---- Extract reads for BLAST  ----
#***********************************

i = 5
genesDetected %>% filter(geneId == genesDetected$geneId[i])


#Filters
minStartLN = 250 #min lenght of start read
minLN = 500 #min read length
nStartPick = 5 #number of reads to pick closest to start (>= minLN)
nDepthRange = 2 #Number of reads to pick closest to start depth (either side, so x2)

blastReads = map_df(genesDetected$geneId, function(myGene){
  #Get a specific GFA file from the masterGFA and simplify it
  myGFA = list(segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId),
               links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId))
  
  myGFA = fixMetacherchant(myGFA)
  simpleGFA = simplifyGFA(myGFA, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
  while(simpleGFA$changed){
    simpleGFA = simplifyGFA(simpleGFA$gfa, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
  }
  
  simpleGFA = simpleGFA$gfa 
  simpleGFA$segments = simpleGFA$segments %>% mutate(depth = KC / LN)
  
  # writeGFA(simpleGFA, "D:/Desktop/mergedGFA.gfa")
  
  #Generate a distance table from reads to start ARG
  distTable = graph_from_data_frame(simpleGFA$links %>% select(from, to))
  distTable = distances(distTable, 
                        v = V(distTable)$name[str_detect(V(distTable)$name, "_start")], 
                        mode = "all")
  
  distTable = distTable %>% apply(2, function(x) ifelse(sum(x == 1) == 2, 0, min(x)))
  distTable = simpleGFA$segments %>% select(-sequence) %>% 
    left_join(data.frame(name = names(distTable), distStart = distTable), by = "name") 
  
  #Select long reads close to the start
  startMatch = distTable %>%
    filter(LN >= minLN, distStart != 0) %>%
    arrange(distStart) %>% 
    mutate(type = "startMatch") %>% 
    slice(1:nStartPick)
  
  #Select long reads with similar depth to start
  depthMatch = distTable %>% filter(distStart == 0 | LN >= minLN) %>%  arrange(depth)
  depthRange = depthMatch %>% filter(distStart == 0) %>% filter(depth == max(depth)) %>% pull(name)
  depthMatch = depthMatch %>% filter(distStart != 0 | name == depthRange[1])
  depthRange = which(depthMatch$name == depthRange[1])
  
  depthMatch = depthMatch[max((depthRange - nDepthRange), 0):min((depthRange + nDepthRange), nrow(depthMatch)),] %>% 
    mutate(type = "depthMatch") %>% filter(distStart != 0)
  
  #Merge all results  and return df
  rbind(depthMatch, startMatch) %>% group_by(name, LN, KC, depth, distStart) %>% 
    summarise(type = ifelse(n() > 1, "bothMatch", type), .groups = 'drop') %>% 
    rbind(distTable %>% filter(distStart == 0, LN >= minStartLN) %>% mutate(type = "ARG")) %>% 
    left_join(simpleGFA$segments %>% select(name, sequence), by = "name") %>% 
    mutate(geneId = myGene)
})

#Generate the FASTA for BLAST
blastReads$blastId = paste0(">", blastReads$geneId, "_", blastReads$name)

write.table(blastReads %>% select(blastId, sequence) %>% pivot_longer(everything()) %>% pull(value), 
            file = paste0(baseFolder, "/temp/", tempName, "/blastReads.fasta"), 
            sep = "\t", quote = F, col.names = F, row.names = F)

# --- run BLAST -----

# ---- Evaluate BLAST output ----
#********************************
blastOutput = read_json(paste0(baseFolder, "/temp/", tempName, "/blastOutput.json"))




