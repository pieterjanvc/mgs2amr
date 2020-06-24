
library(gfaTools)
library(stringr)
library(dplyr)

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

write.csv(genesDetected, paste0(baseFolder, "/temp/", tempName, "/genesDetected.gfa"))


# ---- Extract reads for BLAST  ----
#***********************************
i = 12
genesDetected %>% filter(geneId == genesDetected$geneId[i])

#Get a specific GFA file from the masterGFA and simplify it
myGFA = list(segments = gfa$segments %>% filter(geneId == genesDetected$geneId[i]) %>% select(-geneId),
             links = gfa$links %>% filter(geneId == genesDetected$geneId[i]) %>% select(-geneId))

myGFA = fixMetacherchant(myGFA)
simpleGFA = simplifyGFA(myGFA, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
while(simpleGFA$changed){
  simpleGFA = simplifyGFA(simpleGFA$gfa, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
}

writeGFA(simpleGFA$gfa, "D:/Desktop/mergedGFA.gfa")

#Find the long reads/unitigs near to the start
test = graph_from_data_frame(gfa$links %>% select(from, to))
test = distances(test, 
                 v = V(test)$name[str_detect(V(test)$name, "_start")], 
                 mode = "all")

test = test %>% apply(2, function(x) ifelse(sum(x == 1) == 2, 0, min(x)))


#------------------------


for(myGene in sampleGenes$geneId){ #myGene = sampleGenes$geneId[1]
  #Load the GFA
  gfa = fixMetacherchant(paste0(baseFolder, sample, "/", myGene, ".gfa"))
  
  #Simplify until no more reads to delete or merge
  simpleGFA = simplifyGFA(gfa, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
  while(simpleGFA$changed){
    simpleGFA = simplifyGFA(simpleGFA$gfa, ratioCutOff = 0.3, trimLooseEnds = 60, returnMeta = T, separateStart = T)
  }
  
  #Write new GFA
  writeGFA(simpleGFA$gfa, paste0(baseFolder, sample, "/simplified/", myGene, "_new.gfa"))
  
}

