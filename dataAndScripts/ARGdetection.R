library(stringr)
library(dplyr)
library(purrr)
library(igraph)
library(gfaTools)
library(tidyr)
library(jsonlite)

library(ssh)

# ---- Set parameters ----
#**************************
HOME = T
baseFolder = "./"
tempFolder = ifelse(HOME, "D:/Documents/wslData/meta2amr/temp/EC_AB_1594755817/", "C:/Users/van9wf/Documents/wslData/meta2amrTemp/EC_AB_1594755817/")
#Replace in final version with tempFolder
tempFolderLinux = ifelse(HOME, "/mnt/d/Documents/wslData/meta2amr/temp/EC_AB_1594755817/", "/c/Users/van9wf/Documents/wslData/meta2amrTemp/EC_AB_1594755817/")
tempName = "EC_AB_1594755817"
winToLin = ifelse(HOME, "D:/Documents/wslData/tmp/", "S:/winToLin/")

source(paste0(baseFolder, "dataAndScrips/extraFunctions.R"))

# ---- Merge all MC output ----
#******************************

#cat $(find . -name '*.gfa') > masterGFA.gfa

#Read the master GFA file as gfa object
filePath = paste0(tempFolder, "masterGFA.gfa")
gfa = readGFA(filePath)

#Read the raw file
myFile = str_split(readLines(filePath), "\t")

#Get all IDs for which the file is not empty
message("REMOVE WSL ON LINUX")
baseFolder = "/mnt/d/Documents/wslData/meta2amr"
geneId = system(paste("wsl find ", tempFolder, 
                    " -name '*.gfa' | xargs wc -l"), intern = T)
geneId = str_match(geneId, "(\\d+).*/(\\d+)/graph\\.gfa$")
geneId = data.frame(geneId = geneId[,3], count = as.integer(geneId[,2])) %>% 
  na.omit() %>% filter(count != 0) %>% pull(geneId)


# #The folder number is the gene ID
# geneId = system(paste0("ls -d ",tempFolder, "/*/"), intern = T) %>% 
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
saveRDS(gfa, file = paste0(tempFolder, "masterGFA.rds"))
# gfa = readRDS(paste0(tempFolder, "masterGFA.rds"))


# ---- Detect important ARG ----
#*******************************

#Extract the kmercounts
kmerCounts = gfa$segments %>% 
  mutate(start = str_detect(name, "_start$")) %>% 
  rename(readId = name)

#Get the ARG list
argGenes = read.csv(paste0(baseFolder, "dataAndScripts/argTable.csv"), 
                    colClasses = list(geneId = "character"))

#Detect the most likely genes
genesDetected = kmerCounts %>% filter(start) %>% 
  left_join(argGenes %>% select(geneId, clusterNr, gene, subtype, nBases, name), by = c("geneId" = "geneId")) %>% 
  group_by(geneId, clusterNr, gene, subtype, nBases, name) %>% 
  summarise(length = sum(LN), kmerCount = sum(KC), n = n(), .groups = 'drop') %>% rowwise() %>% 
  mutate(covered = min(1, length / nBases)) %>% 
  group_by(clusterNr) %>% filter(kmerCount == max(kmerCount))

#Make cut-offs based on steepest tangent in curve (maybe not best if very different RA)
cutOff = function(numbers){
  numbers = sort(numbers)
  myData = data.frame(number = numbers[-1],
                      diff = diff(numbers) / diff(1:length(numbers))) 
  myData %>% filter(diff == max(diff)) %>% pull(number) %>% unique()
}

genesDetected = genesDetected %>% filter(covered >= max(0.9, cutOff(genesDetected$covered)))
write.csv(genesDetected, paste0(tempFolder, "genesDetected.csv"), row.names = F)


# ---- Extract reads for BLAST  ----
#***********************************

# genesDetected = read.csv(paste0(tempFolder, "genesDetected.csv"))

blastReads = map_df(genesDetected$geneId, function(myGene){
  cat(paste0("Detecting all GFA paths for ", genesDetected[genesDetected$geneId == myGene, "gene"], " ..."))
  
  #Get the specific GFA file from the masterGFA
  myGFA = list(segments = gfa$segments %>% filter(geneId == myGene) %>% select(-geneId),
               links = gfa$links %>% filter(geneId == myGene) %>% select(-geneId))
  
  myGFA = fixMetacherchant(myGFA)
  # writeGFA(myGFA, "D:/Desktop/mergedGFA.gfa")
  
  #Get the longest start read (reference for now) -- !!! CHANGE LATER
  startRead = myGFA$segments %>% select(name, LN, KC) %>% filter(str_detect(name, "_start$")) %>% 
    filter(LN == max(LN)) %>% pull(name)
  
  pathsTable = pathsToRead(myGFA, startRead) %>%
  mutate(geneId = myGene, id = paste0(">", myGene, "_", startRead, "_", endRead))
  
  cat(" done\n")
  pathsTable
  
})

write.csv(blastReads, paste0(tempFolder, "blastReads.csv"), row.names = F)
# blastReads = read.csv(paste0(tempFolder, "blastReads.csv"))

#Generate a FASTA FILE FROM THE PATHS
fastaPaths = blastReads %>% #filter(geneId == "4967") %>% 
  filter(pathLength > 2500) %>% #??? Should there be cut-off?
  mutate(sequence = str_trunc(sequence, width = 5000, side = "right", ellipsis = "")) %>% # max 5000 bp?
  select(id, sequence) %>% pivot_longer(cols = everything()) %>% select(value)

write.table(fastaPaths, file = paste0(tempFolder, "allPaths.fasta"),
            sep = "\t", quote = F, row.names = F, col.names = F)

#Cluster the reads to reduce BLASTn search
if(HOME){
  system(sprintf("wsl /opt/usearch11.0.667 -cluster_fast %s -id 0.75 -sort size -uc %s",
         paste0(str_replace(tempFolder, "D:/", "/mnt/d/"), "allPaths.fasta"), 
         paste0(str_replace(tempFolder, "D:/", "/mnt/d/"), "allPaths.out")))
} else {
  mySSH = ssh_connect("van9wf@10.200.10.248:22")
  ssh_exec_wait(mySSH, sprintf("/usr/local/usearch/10.0.240/bin/usearch -cluster_fast %s -sort size -id 0.75 -uc %s",
                               paste0(str_replace(tempFolder, "D:/", "/mnt/d/"), "allPaths.fasta"), 
                               paste0(str_replace(tempFolder, "D:/", "/mnt/d/"), "allPaths.out")))
  ssh_disconnect(mySSH)
}

# Use the cluster file to filter the reads to submit to BLAST
fastaPaths = read.table(paste0(str_replace(tempFolder, "D:/", "/mnt/d/"), "allPaths.out"))
blastReads = blastReads %>% filter(pathLength > 2500) %>% 
  mutate(
    sequence = str_trunc(sequence, width = 5000, side = "right", ellipsis = ""),
    id = str_remove(id, ">")
  ) %>% 
  left_join(fastaPaths %>% select(V9, V10), by = c("id" = "V9")) %>% select(-readOrder) %>% 
  distinct() %>% filter(V10 == "*") %>%
  select(id, sequence) %>% mutate(id = paste0(">", id)) %>% pivot_longer(cols = everything()) %>% select(value)

write.table(blastReads, file = paste0(tempFolder, "blastReads.fasta"),
            sep = "\t", quote = F, row.names = F, col.names = F)


# ---- Perform BLASTn search and get results  ----
#*************************************************

#ONLINE SEARCH WITH API
RID = BLASTapi(paste0(tempFolder, "blastReads.fasta"), tempFolder, waitUntilDone = F, secBetweenCheck = 30)

#Unzip the result file (gz)
system(sprintf("wsl cd %s; unzip %s -d %s; cd %s; rm %s; cat * > ../%s; rm -r ../%s", 
               tempFolderLinux, paste0(RID, ".gz"), RID, RID, paste0(RID, ".json"), paste0(RID, ".json"), RID))

# Merge all JSON files into one, and remove the sequences to save space (not needed for now)
myJSON = readLines(paste0(tempFolder, RID, ".json"))
  #Remove seq and other large attributes (not needed)
findLines = str_detect(myJSON, "^\\s+(\"qseq\":)|(\"hseq\":)|(\"midline\":)")
myJSON = myJSON[!findLines]
  #Fix trailing comma on last one before removed
findLines = str_detect(myJSON, "^\\s+\"gaps")
myJSON[findLines] = str_remove(myJSON[findLines], ",$")
  #Update the top names of the JSON files
findLines = str_detect(myJSON, "^\\s+\"BlastOutput2")
myJSON[findLines] = paste0("  \"read", 1:113, "\": {")
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
  #Save
write(myJSON, paste0(tempFolder, RID, "_noSeq.json"))


# ---- Evaluate results and assign bact to ARG  ----
#***************************************************
# RID = "HK1XVM4J016"
# blastOutput = read_json(paste0(tempFolder, RID, "_noSeq.json"))

#Extract the information on the high scoring pairs from the BLAST results
blastOutput = map_df(sapply(myJSON, "[", "report"), function(myResult){
  hits = myResult$results$search$hits
  map_df(hits, function(y){
    cbind(
      map_df(y$hsps, function(x) x[c(1:7, 12)] %>% as_tibble),
      map_df(y$description[1], as_tibble)
    )
  }) %>% mutate(blastId = myResult$results$search$query_title, 
                qLength = myResult$results$search$query_len)
}) %>% mutate(geneId = str_extract(blastId, "^[^_]+"))

test = blastOutput %>% group_by(geneId, blastId, id, accession, title, taxid, sciname, qLength) %>% 
  summarise(bit_score = max(bit_score), evalue = min(evalue), 
            coverage = getCoverage(query_from, query_to, qLength[1]), .groups = "drop") %>% 
  mutate(plasmid = str_detect(title, "plasmid"))


# test2 = test %>% group_by(blastId, plasmid) %>% 
#   filter(coverage >= quantile(test$coverage, 0.9), qLength >= max(2500, quantile(test$qLength, 0.5)))
# 
# test2 = test %>% group_by(blastId, plasmid, sciname) %>% 
#   summarise(n = n(), coverage = max(coverage), readLength = max(qLength), bit_score = max(bit_score))


fastaPaths = fastaPaths %>%
  select(read = V9, cluster = V10) %>%
  mutate(
    cluster = ifelse(cluster == "*", read, cluster),
    readGene = str_extract(read, "^[^_]+"), clusterGene = str_extract(cluster, "^[^_]+"))
fastaPaths = fastaPaths %>% left_join(test, by = c("cluster" = "blastId"))

unique(test$blastId) %in% 
  unique(fastaPaths$cluster)

test2 = data.frame(clust = sort(unique(fastaPaths$cluster)), blast = sort(unique(blastOutput$blastId)))

fastaPaths$bact = str_extract(fastaPaths$sciname, "^[^\\s]+\\s[^\\s]+")

# test3 = fastaPaths %>% filter(readGene == "4672") %>% group_by(read, plasmid) %>% 
#   # filter(coverage >= 0.9) %>% group_by(read) %>% 
#   mutate(score = as.integer(coverage == max(coverage)) + as.integer(bit_score == max(bit_score)),
#          score2 = coverage * bit_score) %>% 
#   filter(score == max(score) | score2 == max(score2)) %>% ungroup() %>% 
#   filter(coverage >= quantile(coverage, 0.90) | bit_score >= quantile(bit_score, 0.90)) %>%
#   mutate(newScore = bit_score / readLength) %>% 
#   group_by(bact) %>% summarise(newScore = sum(newScore)) %>% 
#   filter(newScore >= cutOff(newScore))


test3 = fastaPaths %>% filter(readGene == "4673") %>% group_by(read, plasmid) %>% 
  mutate(score = coverage * bit_score / readLength) %>% 
  group_by(read) %>% mutate(nonDiscr = length(unique(score)) != 1) %>% filter(nonDiscr) %>% group_by(read, plasmid) %>% 
  filter(score == max(score)) %>% ungroup() %>% 
  filter(coverage >= quantile(coverage, 0.90) | score >= quantile(score, 0.90)) %>% distinct() %>% 
  # mutate(newScore = bit_score / readLength) %>% 
  group_by(bact, plasmid) %>%  filter(coverage == max(coverage)) %>% 
  summarise(score = max(score), bit_score = max(bit_score)) %>% arrange(desc(score)) %>% 
  ungroup() %>% filter(score >= cutOff(score))

test3








#---------- OLD

# fastaPaths = read.table(paste0(winToLin, "allPaths.out")) %>% 
#   select(read = V9, cluster = V10) %>% 
#   mutate(
#     cluster = ifelse(cluster == "*", read, cluster),
#     readGene = str_extract(read, "^[^_]+"), clusterGene = str_extract(cluster, "^[^_]+"))
# 
# test = fastaPaths %>% group_by(readGene, clusterGene) %>% summarise(nCommon = n()) %>% 
#   left_join(fastaPaths %>% group_by(readGene) %>% summarise(nReads = n())) %>% 
#   mutate(procent = nCommon / nReads)

#HK1XVM4J016



fastaPath = paste0(tempFolder, "blastReads.fasta")

test = test  %>% 
  left_join(blastReads %>% select(id, sequence) %>% mutate(id = str_remove(id, ">")), by = "id")

test = blastReads %>% group_by(sequence) %>% summarise(n = n())

test2 = map_df(genesDetected$geneId, function(myGene){
  
  fastaPaths = test %>% filter(geneId == myGene) %>% 
    select(id, sequence) %>% pivot_longer(cols = everything()) %>% select(value)
  
  write.table(fastaPaths, file = "D:/Desktop/allPaths.fasta",
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  system("wsl /opt/usearch11.0.667 -cluster_fast /mnt/d/Desktop/allPaths.fasta -id 0.8 -sort size -uc /mnt/d/Desktop/allPaths.out")

  fastaPaths = read.table("D:/Desktop/allPaths.out") %>% filter(V1 == "C")

  test %>% slice(which( str_remove(pathsTable$id, "^>") %in% fastaPaths$V9))  %>%
    select(id, sequence) %>% pivot_longer(cols = everything())
  
})


write.table(test, file = "D:/Desktop/allPaths2.fasta",
            sep = "\t", quote = F, row.names = F, col.names = F)

myGFA = fixMetacherchant("D:/Documents/wslData/meta2amr/temp/EC_AB_1594755817/multi/graph.gfa")
writeGFA(myGFA, "D:/Desktop/mergedGFA.gfa")

baseFolder = "/mnt/d/Documents/wslData/meta2amr"
system(sprintf("wsl %s --tool environment-finder-multi --seq %s --work-dir %s --output %s --env %s",
        "/opt/metacherchant/out/metacherchant.sh",
        paste0(tempFolder, "/multi/checkGene.fasta"),
        paste0(tempFolder, "/multi/"),
        paste0(tempFolder, "/multi/"),
        paste(paste0(tempFolder, "/", genesDetected$geneId, "/graph.txt"), collapse = " ")
        ))

system("wsl /opt/usearch11.0.667 -cluster_fast /mnt/d/Documents/wslData/meta2amr/temp/EC_AB_1594755817/multi/graph.fasta -id 0.75 -uc /mnt/d/Desktop/allPaths.out")

#Generate the FASTA for BLAST
blastReads$blastId = paste0(">", blastReads$geneId, "_", blastReads$name)

write.table(blastReads %>% select(blastId, sequence) %>% pivot_longer(everything()) %>% pull(value), 
            file = paste0(tempFolder, "/blastReads.fasta"), 
            sep = "\t", quote = F, col.names = F, row.names = F)

blastReads$blastId = paste0(blastReads$geneId, "_", blastReads$name)
blastReads$geneId = as.character(blastReads$geneId) 

write.csv(blastReads, paste0(tempFolder, "/blastReads.csv"), row.names = F)

# --- run BLAST -----

# ---- Evaluate BLAST output ----
#********************************
blastOutput = read_json(paste0(tempFolder, "/blastOutput.json"))
RID = "HK1XVM4J016"
blastOutput = read_json(paste0(tempFolder, RID, "_noSeq.json"))

#Extract the information on the high scoring pairs from the BLAST results
blastOutput = map_df(sapply(blastOutput, "[", "report"), function(myResult){
  hits = myResult$results$search$hits
  map_df(hits, function(y){
    cbind(
      map_df(y$hsps, function(x) x[c(1:7, 12)] %>% as_tibble),
      map_df(y$description[1], as_tibble)
    )
  }) %>% mutate(blastId = myResult$results$search$query_title, 
                qLength = myResult$results$search$query_len)
})

#Only keep the best high-scoring pair for each sub-result 
test = blastOutput %>% 
  group_by(blastId) %>%
  filter(evalue == min(evalue)) %>% filter(bit_score == max(bit_score)) %>%
  mutate(bact = str_extract(sciname, "^\\w+\\s\\w+"), plasmid = str_detect(title, "plasmid|Plasmid")) %>% 
  left_join(blastReads %>% select(-sequence), by = "blastId") %>% 
  left_join(argGenes %>% select(geneId, gene), by = "geneId")



#Check the ARG itself
startARG = blastOutput %>% group_by(geneId, gene) %>% filter(type == "ARG") %>% 
  select(geneId, gene, blastId, bit_score, evalue, sciname, bact, plasmid, LN:gene) %>% 
  distinct() %>% group_by(geneId, gene, blastId, plasmid) %>% 
  filter(evalue == max(evalue)) %>% filter(bit_score == max(bit_score))

startARG = startARG %>% group_by(geneId, gene, blastId, sciname, plasmid) %>% 
  summarise(bit_score = max(bit_score), .groups = "drop") %>% 
  group_by(geneId, gene, sciname, plasmid) %>% 
  summarise(n = n(), bit_score = sum(bit_score), .groups = "drop") %>% 
  ungroup() %>% group_by(geneId) %>% filter(n == max(n))


#Check the other reads
otherReads = blastOutput %>% group_by(geneId) %>% filter(type != "ARG") %>% 
  select(geneId, blastId, bit_score, evalue, sciname, bact, plasmid, LN:gene) %>% 
  distinct() %>% group_by(geneId, blastId, plasmid) %>% 
  filter(evalue == max(evalue)) %>% filter(bit_score == max(bit_score))

otherReads = map_df(genesDetected$geneId, function(myGene){
  
  #Next to start 
  nextToStart = otherReads %>% filter(distStart == 1, geneId == myGene) %>% 
    group_by(geneId, gene, bact) %>% summarise(bit_score = sum(bit_score))
 
  if(nrow(nextToStart) == 1){
    otherReads = otherReads %>% filter(distStart == 1, geneId == myGene) %>% 
      group_by(geneId, gene, sciname, plasmid) %>% summarise(n = n(), bit_score = sum(bit_score))
  } else if("bothMatchx" %in% otherReads$type){
    otherReads = otherReads %>% filter(geneId == myGene, type == "bothMatch") %>%
      group_by(geneId, gene, blastId, sciname, plasmid) %>%
      summarise(bit_score = max(bit_score), .groups = "drop") %>%
      group_by(geneId, gene, sciname, plasmid) %>%
      summarise(n = n(), bit_score = sum(bit_score), .groups = "drop") %>%
      ungroup() %>% filter(n == max(n))
  } else {
    otherReads = otherReads %>% filter(geneId == myGene) %>% 
      group_by(geneId, gene, sciname, plasmid) %>%
      summarise(n = n(), bit_score = sum(bit_score), .groups = "drop") %>%
      filter(bit_score == max(bit_score))
  }
  
  return(otherReads)
})

