library(gfaTools)
library(dplyr)
library(tidyr)

#Based off a portion of a gfa (so we know MC will generate output)
gfa = fixMetacherchant("../mixedMetagenomes/EC_KP_MC/1295/graph.gfa")

gfa = readGFA("../mixedMetagenomes/EC_KP_MC/21/graph.gfa")
simple = simplifyGFA(gfa, ratioCutOff = 0.4, trimLooseEnds = 80, separateStart = T)
simple = simplifyGFA(simple, ratioCutOff = 0.4, trimLooseEnds = 80, separateStart = T)

writeGFA(simple, "D:/Desktop/testGraph.gfa")


#Generate fake fastq
genome = readLines("dataAndScripts/testData/sampleStringFastq.txt", 1)

nReads = 250
strLength = nchar(genome)

fastq = lapply(1:nReads, function(x){
  
  start = sample(1:(strLength - 250),1)
  length = sample(c(150,200,250), 1, prob = c(0.1,0.3,0.6))
  
  list(read = substring(genome, start, start + length), 
       qual = paste(sample(c("A", "B", "F", "/", "E", "G", 'H'), length +1, T), collapse = ""))
  
})

#Convert the DF to long format with the values row the lines in the final fastq
fastq = data.frame(
  name = paste0("@read", rep(1:(nReads/2), each = 2), ".", 1:2, " ", rep(1:(nReads/2), each = 2)),
  read = sapply(fastq, "[[", 1), 
  plus = "+",
  qual = sapply(fastq, "[[", 2)) 

fastq = fastq %>% 
  pivot_longer(everything())

#Write the fastq
write.table(fastq$value, "dataAndScripts/testData/fakeFastq.fastq", 
            quote = F, sep = "\t", col.names = F, row.names = F)

