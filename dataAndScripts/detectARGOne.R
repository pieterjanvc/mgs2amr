# detach("package:gfaTools", unload=TRUE)
# library(gfaTools)

baseFolder = "/mnt/meta2amrData/meta2amr/"
database = "/mnt/meta2amrData/pipelineTest/after1200/meta2amr.db"
# database = "/mnt/meta2amrData/pipelineTest/extraTests/extraTests.db"
runId = 0
verbose = 1
pipelineIds = NULL
generateReport = F
forceRedo = F

myId = "6" #452, 845, 18


minBlastLength = 250
#INTEGRATE LATER AS ARGUMENT
Sys.setenv(BLASTDB = "/mnt/meta2amrData/ncbi/blastdb") 
outfmt = "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart send bitscore score length pident nident qcovs qcovhsp"
maxCPU = parallel::detectCores()

#FUNCTIONS
softmax = function(vals, normalise = F, log = T){
  if(normalise) vals = vals / max(vals)
  if(log) vals = log(vals)
  return(exp(vals) / sum(exp(vals)))
}

blast_readOutput = function(file, outfmt, separate = T, includeIssues = F, verbose = 1){
  
  #Load the output csv (zipped)
  blastOut = read.table(file, sep = "\t", quote = "", comment.char = "")
  colnames(blastOut) = strsplit(outfmt, " ")[[1]][-1]
  
  #split multiple matches
  blastOut = blastOut %>% 
    mutate(x = str_count(staxids, ";"), y = str_count(sallacc, ";"), 
           z = x == y | x == 0, plasmid = str_detect(salltitles, "plasmid|Plasmid")) 
  
  if(!includeIssues){
    blastOut = blastOut %>% filter(!(x > 0 & !z))
    issues = data.frame()
  } else {
    #Multiple taxid but not the same number of accessions
    issues = blastOut %>% filter(x > 0 & !z) 
  }
  
  if(nrow(issues) > 0){
    warning(nrow(issues), " rows contain ambiguous accession / taxid results")
  }
  
  if(separate){
    return(bind_rows(
      
      #No merged data
      blastOut %>% filter(y == 0),
      #Identical number of accession and taxids 
      blastOut %>% 
        filter(x > 0 & z) %>% 
        separate_rows(sallacc, staxids, sscinames, sep = ";"),
      #One taxId for multiple accessions
      blastOut %>% 
        filter(x == 0 & y > 0 & z) %>% 
        separate_rows(sallacc, staxids, sscinames, sep = ";"),
      issues
      
    ) %>% select(-x, -y, -z))
  } else {
    return(blastOut %>% select(-x, -y, -z))
  }
  
}


myConn = dbConnect(SQLite(), database)

ARG = dbReadTable(myConn, "ARG") 

#Generate a list out of the settings file
settings = readLines(paste0(baseFolder, "settings.txt"))
settings = settings[str_detect(settings,"^\\s*[^=#]+=.*$")]
settings = str_match(settings, "\\s*([^=\\s]+)\\s*=\\s*(.*)")
settings = setNames(str_trim(settings[,3]), settings[,2])

#Get the AMR prediction models
bactGenomeSize = 
  read.csv(sprintf("%sdataAndScripts/bactGenomeSize.csv",baseFolder))

toProcess = dbReadTable(myConn, "pipeline") %>% 
  filter(statusCode > 2)

sampleIndex = which(toProcess$pipelineId == myId)
dbDisconnect(myConn)

sample = toProcess$tempFolder[sampleIndex]
myPipelineId = toProcess$pipelineId[sampleIndex]
sampleName = str_extract(sample, "[^\\/]+(?=_\\d+$)")

inputfileBP = grep("sequences length",
                   readLines(sprintf("%s/metacherchant_logs/log", sample)),
                   value = T, fixed = T) %>%
  str_extract("([\\d'])+(?=\\s\\()") %>%
  str_replace_all("'", "") %>% as.numeric()

#Grab the detected ARG from the previous step
myConn = dbConnect(SQLite(), database,synchronous = NULL)
sqliteSetBusyHandler(myConn, 30000)

genesDetected = tbl(myConn, "detectedARG") %>% 
  filter(pipelineId == myPipelineId) %>% as.data.frame() %>% 
  mutate(cover = ifelse(type == 'noFragments', cover1, cover2))
dbDisconnect(myConn)

newLogs = data.frame(
  timeStamp = as.integer(Sys.time()), 
  actionId = 1, actionName = "Start Annotation & Prediction")

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    sprintf("Start Annotation & Prediction for pipelineId %i (%i/%i)\n", 
                            myPipelineId, sampleIndex, nrow(toProcess)))}

sample

# ---- LOAD DATA ----
#*************************
if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Loading BLASTn output ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 2, 
                              "Start loading BLASTn output"))

#Read all blast output
blastOut = lapply(
  list.files(sample, full.names = T,
             pattern = "blastSegmentsClustered\\d+.csv.gz"),
  blast_readOutput, outfmt = outfmt) %>% bind_rows() %>% 
  mutate(geneId = str_extract(qseqid, "^([^_]+)"))


# # ---- RERUN BLAST FOR IDENTICAl RESULTS ----
# #********************************************
# if(!(file.exists(sprintf("%s/expand_1.csv.gz", sample)) | 
#    file.exists(sprintf("%s/expand_2.csv.gz", sample))) | forceRedo){
#   
#   #Get segments that have identical blast results (thus might miss some because of 250 lim)
#   rerun = blastOut %>% 
#     select(segmentId = qseqid, bitscore, geneId, nident, 
#            qlen, length) %>% 
#     group_by(segmentId, geneId) %>% 
#     summarise(x = n_distinct(bitscore) == 1, .groups = "drop",
#               identity = nident[1] / qlen[1],
#               cover = length[1] / qlen[1]) %>% 
#     filter(x)
#   
#   #Set these general blast args
#   blastArgs = list(
#     db = "nt",
#     evalue = "1e-20",
#     word_size = 64,
#     max_target_seqs = 5000, #Set max of 5000 results per target
#     max_hsps = 1,
#     qcov_hsp_perc = 100, #Only consider perfect matches
#     perc_identity = 100, #Only consider perfect matches
#     taxidlist = sprintf("%s/bact.txids", sample) #limit to taxa found in first part
#   )
#   
#   write(unique(blastOut$taxid), sprintf("%s/bact.txids", sample), sep = "\n")
#   
#   #Blast all perfect matches together
#   toFasta = rerun %>% filter(cover == 1, identity == 1)
#   
#   if(nrow(toFasta) > 0){
#     
#     myFasta = map_df(1:length(list.files(sample, pattern = "blastSegmentsClustered\\d.csv.gz")), function(x){
#       fasta_read(sprintf("%s/blastSegmentsClustered%i.fasta", sample, x),
#                  type = "n") %>% filter(id %in% toFasta$segmentId)})
#     
#     fasta_write(myFasta$seq, sprintf("%s/expand_1.fasta", sample), 
#                 myFasta$id, type = "n")
#     
#     #Run local blastn on the new database
#     system(sprintf('blastn -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -num_threads %i -qcov_hsp_perc %.2f -perc_identity %.2f -taxidlist %s -outfmt "%s" | gzip > "%s"',
#                    blastArgs$db, sprintf("%s/expand_1.fasta", sample),
#                    blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
#                    blastArgs$max_hsps, maxCPU, blastArgs$qcov_hsp_perc,
#                    blastArgs$perc_identity, blastArgs$taxidlist, outfmt,
#                    sprintf("%s/expand_1.csv.gz", sample)))
#   }
#   
#   #Imperfect ones get blased separately
#   toFasta = rerun %>% filter(cover < 1 | identity < 1)
#   if(nrow(toFasta) > 0){
#     
#     myFasta = map_df(1:length(list.files(sample, pattern = "blastSegmentsClustered\\d+.csv.gz")), function(x){
#       fasta_read(sprintf("%s/blastSegmentsClustered%i.fasta", sample, x),
#                  type = "n") %>% filter(id %in% toFasta$segmentId)})
#     
#     fasta_write(myFasta$seq, sprintf("%s/expand_2.fasta", sample), 
#                 myFasta$id, type = "n")
#     
#     #Set the settings to the min cover and identity of imperfect ones
#     blastArgs$qcov_hsp_perc = floor(min(toFasta$cover) * 10000) / 10000
#     blastArgs$perc_identity = floor(min(toFasta$identity) * 10000) / 10000
#     
#     #Run local blastn on the new database
#     system(sprintf('blastn -db "%s" -query "%s" -task megablast -evalue %s -word_size %i -max_target_seqs %i -max_hsps %i -num_threads %i -qcov_hsp_perc %.2f -perc_identity %.2f -taxidlist %s -outfmt "%s" | gzip > "%s"',
#                    blastArgs$db, sprintf("%s/expand_2.fasta", sample),
#                    blastArgs$evalue, blastArgs$word_size, blastArgs$max_target_seqs, 
#                    blastArgs$max_hsps, maxCPU, blastArgs$qcov_hsp_perc,
#                    blastArgs$perc_identity, blastArgs$taxidlist, outfmt,
#                    sprintf("%s/expand_2.csv.gz", sample)))
#   }
#   
#   
# }

expandBlast = list.files(sample, full.names = T, pattern = "expand_\\d+.csv.gz")

if(length(expandBlast) > 0){
  
  #Get the new results
  blastOut2 = lapply(
    list.files(sample, full.names = T,
               pattern = "expand_\\d+.csv.gz"),
    blast_readOutput, outfmt = outfmt) %>% bind_rows() %>% 
    mutate(geneId = str_extract(qseqid, "^([^_]+)"))
  
  blastOut = bind_rows(
    blastOut %>% filter(!qseqid %in% unique(blastOut2$qseqid)),
    blastOut2
  )
  
  rm(blastOut2)
}

#FIX FROM HERE
# backup1 = blastOut 
# blastOut = backup1

#Extract data we need + transform
blastOut = blastOut %>% 
  select(query_title = qseqid, taxid = staxids, accession = sallacc, 
         bact = sscinames, plasmid, bit_score = bitscore, score, 
         identity = nident, query_len = qlen, query_from = qstart, query_to = qend, 
         hit_from = sstart, hit_to = send, align_len = length, geneId) %>% 
  mutate(bact = str_remove_all(bact, "[^\\w\\s]"),
         geneId = as.integer(geneId), taxid = as.integer(taxid)) %>% 
  extract(bact, c("genus", "species", "extra"), 
          regex = "(\\w+)\\s+(\\w+)($|\\s+.*)") %>% 
  #Sp. will be pasted with taxid to make it unique
  mutate(species = ifelse(species == "sp", paste0(species, taxid), species)) %>% 
  filter(!genus %in% c("uncultured", "Uncultured", "mixed", "Bacterium"),
         !species %in% c("bacterium", "Bacterium") & 
           !is.na(species)) %>% 
  #Extract subspecies, strain and plasmid info
  mutate(
    subspecies = str_extract(extra, "(?<=subsp\\s)[^\\s]+"),
    strain = str_extract(extra, "(?<=strain\\s)[^\\s]+") %>% 
      str_remove(" chromosome.*| plasmid.*| complete.*")
  )

#Expand results from clustering segments before blast
clusterOut = read.table(paste0(sample, "/blastSegments.out")) %>% 
  select(segmentId = V9, clusterId = V10) %>% distinct() %>% 
  mutate(
    clusterId = ifelse(clusterId == "*", segmentId, clusterId)
  ) %>% 
  filter(clusterId %in% blastOut$query_title) 

blastOut = clusterOut %>% 
  left_join(blastOut, by = c("clusterId" = "query_title")) %>% 
  #Get the coverage and identity of each alignment
  rowwise() %>% 
  mutate(
    ident = identity / align_len,
    coverage = min(align_len / query_len, 1),
  ) %>% ungroup() %>% 
  mutate(geneId = str_extract(segmentId, "^([^_]+)"))

#Fix some known issues with naming
blastOut[blastOut$genus == "Enterobacter" & 
           blastOut$species == "aerogenes","genus"] = "Klebsiella"


#--
# backup2 = blastOut
# blastOut = backup2
#--

if(verbose > 0){cat("done\n")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 3, 
                              "Finished loading BLASTn output"))


# ---- MAPPING DATA TO GFA ----
#******************************

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Mapping alignments to GFA ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 4, 
                              "Start mapping"))


segmentsOfInterest = read.csv(paste0(sample, "/segmentsOfInterest.csv"))

#Get information about the blasted seq where they are in the path of GFA
pathData = list.files(
  paste0(sample, "/genesDetected/simplifiedGFA"), 
  pattern = ".gfa", full.names = T)

pathData = map_df(pathData, function(myGFA){
  
  #Read the GFA
  fullGFA = gfa_read(myGFA)
  
  myGene = str_extract(fullGFA$segments$name[1], "^\\d+")
  
  segmentOfInterest = segmentsOfInterest %>% 
    filter(geneId == myGene, type != "fragmentsOnly") 
  
  map_df(1:nrow(segmentOfInterest), function(i){
    
    myGroup = segmentOfInterest %>% filter(name == segmentOfInterest$name[i]) %>% 
      pull(group)
    
    gfa = gfa_filterSegments(fullGFA, segments = (fullGFA$segments %>% 
                               filter(group == myGroup) %>% pull(name)))
    
    0#If there are linked segments ...
    if(nrow(gfa$links) > 0){
      
      # #Get the (largest) start segment
      # segmentOfInterest = gfa$segments %>% 
      #   filter(str_detect(name, "_start$")) %>%
      #   filter(LN == max(LN)) %>% filter(KC == max(KC)) %>% 
      #   dplyr::slice(1) %>% pull(name)
      # 
      # if(length(segmentOfInterest) == 0){
      #   return(data.frame())
      # }
      
      #Get all semenents in path to start
      x = gfa_pathsToSegment(gfa, segmentOfInterest$name[i], returnList = T, 
                             pathSegmentsOnly = T, verbose = F) %>% 
        map_df(function(path){
          data.frame(
            pathId = path$id,
            orientation = path$orientation,
            pathType = path$pathType,
            segmentId = path$segmentOrder,
            dist = path$dist,
            LN = path$LN,
            KC = path$KC)
        }) %>% 
        mutate(geneId = str_extract(segmentId, "^\\d+"))
      
      
      #If no paths, return empty frame
      if(nrow(x) == 0 ){
        return(data.frame())
      }
      #Otherwise add the path data 
      else{
        
        
        x = x %>% 
          # left_join(
          #   pathsToSegmentTable(gfa, segmentOfInterest) %>% 
          #     filter(dist < Inf) %>% group_by(segment) %>% 
          #     filter(dist == min(dist)) %>% 
          #     select(segment, segmentEnd, KC, LN, dist) %>% 
          #     mutate(geneId = str_extract(segment, "^\\d+")) %>% 
          #     distinct(),
          #   by = c("segmentId" = "segment")
          # ) %>% 
          # filter(LN >=250) %>%
          group_by(geneId, pathId) %>% mutate(
            pathId = as.integer(pathId),
            order = n():1,
            #These files contain full GFA structures, no fragments
            type = segmentOfInterest$type[i]) %>% 
          ungroup()
        
        # if(any(x$pathType == 2)){
        #   
        #   #In case of circular loop, only keep 1 path that has least competition
        #   toRemove = x %>% filter(pathType != 2) %>% group_by(pathId, orientation) %>% 
        #     summarise(LN = sum(LN) - 30*(n() - 1), .groups = "drop") %>% 
        #     filter(LN == max(LN)) %>% dplyr::slice(1) %>% pull(orientation)
        #   
        #   x = x %>% filter(!(pathType == 2 & orientation == toRemove))
        # }
        
        return(x)
        
      }
      
      
      
    } else if(nrow(gfa$segments) > 0){
      #Case where there are no links, just one or more segments
      gfa$segments %>% 
        select(segmentId = name, KC, LN) %>% 
        mutate(
          pathId = 1:n(),
          orientation = 0,
          geneId = str_extract(segmentId, "^\\d+"),
          dist = ifelse(str_detect(segmentId, "_start$"), -LN, 0),
          order = 1,
          type = "full"
        )
      
    } else{
      #File is empty
      return(data.frame())
    }
  })
  
})


#Calculate the score to bit_score conversion factor for each gene
extraBits = blastOut %>% filter(ident == 1, coverage == 1) %>% 
  group_by(geneId) %>% 
  summarise(bitConst = mean(bit_score / score), .groups = "drop")

extraBits = pathData %>% 
  filter((LN < minBlastLength & dist >= 0 ) | (dist < 0 & LN < 75)) %>% 
  group_by(geneId, pathId) %>% 
  summarise(score = sum(LN) - 30*n(), .groups = "drop") %>% 
  left_join(extraBits, by = "geneId") %>% 
  mutate(bitScore = bitConst * score)

#Only work with segments >= 250 in paths
pathData = pathData %>% 
  filter(LN >= minBlastLength | (dist < 0 & LN > 74)) %>% 
  group_by(geneId, pathId) %>% 
  arrange(geneId, pathId, desc(order)) %>% 
  mutate(order = n():1) %>% 
  ungroup()

#Also add the fragment data to the paths
fragments = gfa_read(paste0(sample, "/fragmentGFA.gfa"))

if(nrow(fragments$segments) > 0){
  #Fragments were merged, so treated as a non-start piece
  fragments = fragments$segments %>% 
    select(segmentId = name, LN, KC, geneId) %>% 
    group_by(geneId) %>% 
    mutate(pathId = 1, orientation = 0:(n()-1),
           order = 1, type = "fragment", dist = -1)
} else {
  fragments = data.frame()
}

pathData = bind_rows(pathData, fragments)  %>% group_by(geneId) %>% 
  mutate(orientation = ifelse(type == "fragment" & any(type == "full"), 
                              abs(orientation[type == "full"][1]-1), orientation)) %>% 
  ungroup()


pathData = pathData %>% group_by(geneId, pathId, orientation) %>% 
  mutate(
    depth = KC / LN,
    totalLN = sum(LN[dist >= 0]) - 30*sum(dist >= 1)) %>% ungroup()

# #Only keeps paths that have results from blast
# pathData = pathData %>% left_join(data.frame(
#   segmentId = unique(blastOut$segmentId),
#   result = T
# ),by = "segmentId") %>% 
#   group_by(geneId, pathId) %>% 
#   filter(!any(is.na(result))) %>% ungroup() %>% select(-result)


if(verbose > 0){cat("done\n")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 5, 
                              "Finished mapping"))


#Add the grouping info to detected genes
genesDetected = genesDetected %>%
  left_join(ARG %>% select(geneId, gene, subtype), by = "geneId") %>%
  select(pipelineId, geneId, gene, subtype, everything()) %>% 
  mutate(ARGgroup = 1:n())

# ---- Filter / group bacteria ----
#**********************************

if(verbose > 0){cat(format(Sys.time(), "%H:%M:%S -"), 
                    "Group overlapping bacterial species ... ")}
newLogs = rbind(newLogs, list(as.integer(Sys.time()), 8, 
                              "Start grouping bacteria"))

blastSegments = 
  read.csv(paste0(sample, "/blastSegments.csv")) %>% 
  select(name, LN, KC) %>% 
  mutate(start = ifelse(str_detect(name, "_start$"), T, F))

blastOut = blastOut %>% 
  left_join(blastSegments, by = c("segmentId" = "name")) %>% 
  mutate(KC = KC * coverage * ident)


# #fix missing start wiht fragments (change later in blastPrep naming !!!)
# blastOut$start[blastOut$segmentId %in% fragments$segmentId] = T
# blastOut$segmentId[blastOut$segmentId %in% fragments$segmentId] = 
#   paste0(blastOut$segmentId[blastOut$segmentId %in% fragments$segmentId], "_start")
# pathData$segmentId[pathData$segmentId %in% fragments$segmentId] = 
#   paste0(pathData$segmentId[pathData$segmentId %in% fragments$segmentId], "_start")


#Make sure that taxid always has the same genus / species name and vice versa
blastOut = blastOut %>% 
  mutate(bact = paste(genus, species)) %>%
  group_by(taxid, genus, species) %>% 
  mutate(myCount = n()) %>% 
  group_by(taxid) %>% 
  mutate(genus = genus[myCount == max(myCount)][1],
         species = species[myCount == max(myCount)][1]) %>% 
  group_by(taxid, genus, species) %>% 
  mutate(myCount = n()) %>% 
  group_by(genus, species) %>% 
  mutate(taxid = taxid[myCount == max(myCount)][1]) %>% 
  ungroup()

# backup = blastOut
# blastOut = backup
# if(nrow(fragments) > 0){
#   fragments = fragments %>% 
#     mutate(geneId = as.integer(geneId)) %>% 
#     left_join(genesDetected %>% select(geneId, fragType = type), by = "geneId")
# }




blastOut$start[blastOut$start & 
                 ! blastOut$segmentId %in% segmentsOfInterest$from] = F

#Only consider segmenst in the graph that are also the same distance
#in the actual geome alignments
tempVar = blastOut %>% 
  group_by(geneId, accession) %>% 
  filter(any(start)) %>% 
  mutate(sMin = ifelse(hit_from < hit_to, hit_from, hit_to),
         sMax = ifelse(hit_from >= hit_to, hit_from, hit_to),
         after = ifelse(sMin > sMin[start], T, F),
         dist = ifelse(after, sMin - sMax[start], sMin[start] - sMax) + 29) %>%
  ungroup() %>% 
  left_join(pathData %>% select(segmentId, dist1 = dist) %>% 
              distinct(), by = "segmentId") %>% 
  rowwise() %>% 
  #Repeats can make segments shorter because loop is cut off, so provide some slack
  mutate(correct = between(dist, dist1 -500, dist1 + 500) | start) %>% 
  filter(correct)

wrongLocation = blastOut %>% filter(!segmentId %in% tempVar$segmentId)
# tempVar = blastOut %>% filter(!segmentId %in% wrongLocation$segmentId)
blastOut = bind_rows(
  tempVar,
  blastOut %>% 
    filter(geneId %in% (genesDetected %>% 
                          filter(!geneId %in% segmentsOfInterest$geneId) %>% 
                          pull(geneId)))
  )
rm(tempVar)

allBact = blastOut %>%
  # filter(geneId %in% c("4793")) %>%
  select(segmentId, geneId, bit_score, coverage,
         accession, taxid, genus, species, extra, plasmid, KC, LN, start) %>%
  group_by(segmentId, geneId, accession) %>%
  filter(bit_score == max(bit_score)) %>% 
  dplyr::slice(1) %>% ungroup() %>%
  mutate(pathScore = bit_score * coverage) %>%
  # mutate(pathScore = bit_score) %>% 
  group_by(segmentId, accession, plasmid) %>% 
  filter(pathScore == max(pathScore)) %>% dplyr::slice(1) %>% ungroup() %>% 
  left_join(pathData %>% 
              select(pathId, segmentId, order, orientation) %>% 
              mutate(orientation = ifelse(str_detect(segmentId, "_start$"),
                                          -1, orientation)) %>% 
              distinct(), by = "segmentId") %>% 
  filter(!is.na(order))  %>% 
  arrange(geneId, pathId, accession, plasmid, order) %>% 
  group_by(geneId, pathId, accession, plasmid) %>% 
  filter(order == 1:n()) %>% 
  
  
  # group_by(geneId, pathId, accession, plasmid, orientation) %>% 
  # filter(all(2:max(order) %in% order) | all(order == 1)) %>% 
  
  
  # group_by(segmentId, accession) %>% filter(pathScore == max(pathScore)) %>% 
  # dplyr::slice(1) %>% group_by(geneId, accession, pathId) %>% 
  # mutate(fullPath = sum(pathScore))  
  
  #Calculate path score
  group_by(geneId, accession, plasmid, orientation, pathId) %>%
  mutate(fullPath = sum(pathScore)) %>%
  #Pick best path for each orientation
  group_by(geneId, accession, plasmid, orientation) %>%
  filter(pathId == pathId[fullPath == max(fullPath)][1]) %>%
  #Sum the 3 orientations to get full path score
  group_by(geneId, accession, plasmid) %>%
  mutate(fullPath = sum(pathScore)) %>%
  group_by(geneId, accession) %>%
  filter(plasmid == plasmid[fullPath == max(fullPath)][1]) %>% 
  ungroup()

#Calculate the exta bits to add to the total path score
test = allBact %>% filter(orientation != -1) %>% 
  select(geneId, pathId, accession) %>% distinct() %>% 
  left_join(extraBits, by = c("geneId", "pathId")) %>%
  group_by(geneId, accession) %>% 
  mutate(score = sum(score, na.rm = T),
         bitScore = sum(bitScore, na.rm = T)) %>% 
  group_by(geneId, accession) %>% 
  summarise(score = max(score), bitConst = max(0, bitConst, na.rm = T),
                        bitScore = max(bitScore), .groups = "drop")
#Add the extra bits
allBact = allBact %>% 
  left_join(test, by = c("geneId", "accession")) %>%
  mutate(score = replace_na(score, 0), bitScore = replace_na(bitScore, 0)) %>% 
  mutate(fullPath = fullPath + bitScore)
  
allBact = allBact %>% 
  left_join(
    allBact %>% group_by(geneId, accession) %>% 
      summarise(fullPath = max(fullPath), .groups = "drop") %>% 
      group_by(accession) %>% 
      summarise(fullPath = sum(fullPath),.groups = "drop") %>% 
      transmute(accession, rank = rank(fullPath)),
    by = "accession"
  ) %>% 
  group_by(geneId, accession) %>% 
  mutate(path0 = pathId[orientation == 0][1],
         maxOrder0 = max(0,order[orientation == 0]),
         path1 = pathId[orientation == 1][1],
         maxOrder1 = max(0,order[orientation == 1]),
         LN = sum(LN), KC = sum(KC)) %>% 
  group_by(geneId, taxid, pathId) %>% 
  filter(fullPath == max(fullPath)) %>%
  filter(accession %in% accession[rank == max(rank)][1]) %>%
  # filter(rank == max(rank)) %>% dplyr::slice(1) %>% 
  ungroup() 


# test = allBact %>% group_by(geneId, accession) %>% 
#   summarise(fullPath = max(fullPath), .groups = "drop") %>% 
#   group_by(accession) %>% 
#   summarise(fullPath = sum(fullPath),.groups = "drop") %>% 
#   transmute(accession, rank = rank(fullPath))

#Remove genes with no blast data
# geneIds = genesDetected %>% pull(geneId)
# geneIds = geneIds[geneIds %in% unique(allBact$geneId)]

## Check the overlap between segments of different species
## if one is completely contained within another (and smaller), it's a duplicate
# allBact = allBact %>% group_by(geneId, path0, path1) %>% 
#   filter(maxOrder0 == max(maxOrder0), maxOrder1 == max(maxOrder1)) %>%
#   group_by(geneId, path1) %>% filter(!(is.na(path0) & max(maxOrder0) > 0)) %>% 
#   group_by(geneId, path0) %>% filter(!(is.na(path1) & max(maxOrder1) > 0))

allBact = allBact %>% 
  group_by(geneId, accession, taxid, genus, species, plasmid, LN, KC) %>% 
  summarise(extension = sum(pathScore[!start]) + bitScore[1],
            across(c(fullPath:maxOrder1), max), .groups = "drop") %>% 
  mutate(
    depth = KC / LN,
    geneId = as.integer(geneId)) %>% 
  left_join(
    genesDetected %>%
      select(geneId, gene, subtype, cover, type), 
    by = "geneId"
  ) %>% distinct() %>% 
  mutate(pathScore = fullPath)

# #Check the overlap between segments of different species
# # if one is completely contained within another (and smaller), it's a duplicate
# allBact = map_df(geneIds, function(geneId){
#   
#   #Examine the paths of one ARG
#   myGene = allBact %>% filter(geneId %in% {{geneId}})
# 
#   #Get all the segments per path
#   bactSegments = myGene %>% 
#     group_by(segmentId, taxid, pathId) %>% 
#     summarise(.groups = "drop") %>% 
#     group_by(taxid, pathId) %>% 
#     summarise(x = list(segmentId), .groups = "drop") 
#   
#   if(length(bactSegments$taxid) == 0){
#     
#     #There is no blast output for this gene
#     return(data.frame())
#     
#   } else if(length(bactSegments$taxid) == 1){
#     
#     #There is only one taxid for the gene
#     myGene = myGene %>% 
#       group_by(geneId, taxid, genus, species, plasmid) %>% 
#       summarise(n = n(), 
#                 extension = sum(pathScore[!start]) + bitScore[1],
#                 pathScore = sum(pathScore) + bitScore[1], 
#                 KC = sum(KC), LN = sum(LN),
#                 .groups = "drop") %>%
#       group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
#       arrange(desc(pathScore))
#     
#     myTaxid = bactSegments$taxid
#     
#     return(myGene %>% filter(taxid %in% myTaxid) %>% 
#              group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
#              dplyr::slice(1) %>% ungroup())
#     
#   } else {
#     
#     #There are multiple taxids for the gene
#     #Check per path if there are redundant taxids (contained within other)
#     pathResult = map_df(unique(bactSegments$pathId), function(pathId){
#       
#       #Get segments for the path
#       pathSegments = bactSegments %>% filter(pathId == {{pathId}})
#       
#       if(nrow(pathSegments) > 1){
#         #Build a matrix for species overlap per segment
#         gfaMatrix = sapply(pathSegments$x, function(x){
#           len = length(x)
#           sapply(pathSegments$x, function(y){
#             sum(x %in% y) / len
#           })
#         })
#         
#         rownames(gfaMatrix) = pathSegments$taxid
#         colnames(gfaMatrix) = pathSegments$taxid
#       }
#       
#       myPath = myGene %>% 
#         filter(taxid %in% pathSegments$taxid) %>% 
#         group_by(geneId, taxid, genus, species, plasmid) %>% 
#         # summarise(n = n(), 
#         #           extension = sum(pathScore[!start]),
#         #           pathScore = sum(pathScore), 
#         #           KC = sum(KC), LN = sum(LN),
#         #           .groups = "drop") %>%
#         summarise(n = n(), 
#                   extension = sum(pathScore[!start]) + bitScore[1],
#                   pathScore = sum(pathScore) + bitScore[1], 
#                   KC = sum(KC), LN = sum(LN),
#                   .groups = "drop") %>% 
#         group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
#         arrange(desc(pathScore))
#       
#       myTaxid = myPath$taxid %>% unique()
#       i = 1
# 
#       #Remove taxid that are contained within others
#       while(i < length(myTaxid)){
#         myIndex = which(colnames(gfaMatrix) == myTaxid[i])
#         duplicates = data.frame(
#           name = colnames(gfaMatrix),
#           row = gfaMatrix[myIndex,],
#           col = gfaMatrix[,myIndex]
#         ) %>% filter(row == 1, col != 1) %>% pull(name)
#         
#         myTaxid = setdiff(myTaxid, duplicates)
#         i = i + 1
#       }
#       
#       return(myPath %>% filter(taxid %in% myTaxid) %>% 
#                group_by(taxid) %>% filter(pathScore == max(pathScore)) %>% 
#                dplyr::slice(1) %>% ungroup())
# 
#     })
#     
#     return(pathResult)
#     
#   } 
#   
# }) %>%  mutate(
#     depth = KC / LN,
#     geneId = as.integer(geneId)) %>% 
#   left_join(
#     genesDetected %>%
#       select(geneId, gene, subtype, cover, type), 
#     by = "geneId"
#   ) %>% distinct()

#Adjust the scored for genes by presences of other genes in the same bact
# myData = genomeARG
# myData$taxid = myData$accession

adjustBact = function(myData, bactGroupStart = 0){
  
  if(nrow(myData) == 0) {
    return(list(myClusters = data.frame(), bactList = data.frame()))
  }
  
  if(n_distinct(myData$geneId) == 1 | n_distinct(myData$taxid) == 1){
    if(n_distinct(myData$taxid) == 1){
      
      #There is only one taxid
      bactList = myData %>% 
        select(taxid, genus, species, depth, pathScore) %>% 
        mutate(bactGroup = 1 + bactGroupStart, prob = 1, val = NA,
               runId, pipelineId = myPipelineId) %>% 
        filter(depth == max(depth)) %>% dplyr::slice(1) %>% 
        left_join(bactGenomeSize %>% select(genus, species, size) %>% distinct(),
                  by = c("genus", 'species')) 
      
      return(list(
        myClusters = data.frame(
          bactGroup = 1 + bactGroupStart,
          geneId = myData$geneId[1]),
        bactList = bactList))
    } else {
      
      #There is only one gene but multiple taxid
      geneMatrix = matrix(myData$pathScore / max(myData$pathScore), ncol = 1)
      colnames(geneMatrix) = unique(myData$geneId)
      rownames(geneMatrix) = myData$taxid
      myClusters = data.frame(
        taxid =  myData$taxid,
        bactGroup = 1 + bactGroupStart,
        primary = geneMatrix[,1] == 1,
        geneId = as.integer(unique(myData$geneId)),
        val = geneMatrix[,1])
    }
    
    
  } else {
    
    #There are multiple genes
    geneMatrix = myData %>%
      mutate(geneId = as.integer(geneId)) %>% 
      group_by(geneId, bact = taxid) %>%
      summarise(pathScore = round(max(pathScore),2), 
                extension = round(max(extension),2), 
                .groups = "drop") %>%
      left_join(genesDetected %>%
                  select(geneId, ARGgroup), by = "geneId")
    
    geneMatrix = geneMatrix %>%
      select(bact, geneId, pathScore) %>%
      pivot_wider(bact,
                  names_from = "geneId", values_from = "pathScore",
                  values_fill = 0) %>%
      column_to_rownames("bact") %>% as.matrix()
    
    #Normalise the matrix per column
    geneMatrix = apply(geneMatrix, 2, function(x) x / max(x))
    
    # Adjust by cover
    myAdjustment = data.frame(
      geneId = as.integer(colnames(geneMatrix))
    ) %>% left_join(
      genesDetected %>% select(geneId, cover), by = "geneId")
    
    #-----------
    # During adjustment make sure to weigh the influence of other ARG based
    # on the similarity in coverage. If the coverage is very different, the
    # weight of adjustment in lower (prevents FP association of low cover to high)
    ncols = ncol(geneMatrix)
    myAdjustment = (1 - abs(matrix(myAdjustment$cover, nrow = ncols, ncol = ncols) -
                              matrix(myAdjustment$cover, nrow = ncols, ncol = ncols, byrow = T)))
    
    
    myAdjustment = sapply(1:ncols, function(col){
      t(t(geneMatrix[,-col]) * myAdjustment[col,-col]) %>% rowSums() / 
        (ncols - 1)
    }) + 1
    
    #---------------
    ##Old adjustment techinique
    # geneMatrix = t(t(geneMatrix) * myAdjustment$cover)
    # 
    # #Calulcate the adjustment for each cell based on other scores in the row
    # # the higher the scores (e.g. other gened detected), the more scaled up
    # myAdjustment = (matrix(rowSums(geneMatrix), nrow = nrow(geneMatrix), 
    #                        ncol = ncol(geneMatrix)) - geneMatrix) /
    #   (ncol(geneMatrix) - 1) + 1
    
    #----------------
    
    geneMatrix = geneMatrix * myAdjustment
    geneMatrix = apply(geneMatrix, 2, function(x) x / max(x))
    
    
    myClusters = as_tibble(geneMatrix > 0, rownames = "taxid") %>%
      group_by(across(c(-taxid))) %>% 
      mutate(bactGroup = cur_group_id()) %>% ungroup()
    
    myClusters = as_tibble(geneMatrix, rownames = "taxid") %>% 
      left_join(myClusters %>% select(taxid, bactGroup), by = "taxid") 
    myClusters$primary = apply(myClusters[,!colnames(myClusters) %in% c("taxid", "bactGroup")], 1,
                               function(x) any(x == 1))
    
    #new ---
    x = myClusters %>% group_by(bactGroup) %>% filter(!all(!primary)) %>% 
      ungroup() 
    
    
    y = myClusters %>% group_by(bactGroup) %>% filter(all(!primary)) %>% 
      ungroup()
    # mutate(across(c(-taxid, -bactGroup, -primary), function(x) x > 0))
    
    
    a = apply(x %>% select(-taxid, -bactGroup, -primary), 1, function(z){
      
      apply(y %>% select(-taxid, -bactGroup, -primary), 1, function(a) cor(z,a))
      
    })
    
    colnames(a) = x$bactGroup
    
    y = cbind(a,y)
    y = y[,1:which(colnames(y) == "taxid")]
    
    y = y %>% pivot_longer(-taxid, names_to = "cluster") %>% 
      group_by(taxid) %>% filter(value == max(value)) %>% ungroup() %>% 
      mutate(cluster = str_remove(cluster, "\\.\\d+$"))
    
    myClusters = myClusters %>% left_join(y, by = "taxid") %>% 
      mutate(bactGroup = ifelse(is.na(cluster), bactGroup, cluster)) %>% 
      select(-cluster)
    #-----
    
    myClusters = myClusters %>% group_by(bactGroup) %>% filter(!all(!primary)) %>% 
      ungroup() %>% mutate(bactGroup = bactGroupStart + as.factor(bactGroup) %>% as.integer()) %>% 
      pivot_longer(c(-taxid, -bactGroup, -primary, -value), names_to = "geneId", values_to = "val") %>% 
      mutate(geneId = as.integer(geneId), taxid = as.integer(taxid)) %>% 
      filter(val > 0 | primary)
    
  }
  
  bactList = myClusters %>% group_by(bactGroup, taxid) %>% 
    summarise(prob = sum(val), .groups = "drop") %>% 
    group_by(bactGroup) %>% 
    mutate(prob = softmax(prob)) %>% ungroup() %>%
    mutate(taxid = as.integer(taxid), bactGroup = as.integer(bactGroup)) %>% 
    left_join(allBact %>% select(taxid, genus, species) %>% 
                distinct(), by = "taxid") %>% 
    group_by(bactGroup) %>% 
    filter(prob >= min(
      sort(unique(prob), decreasing = T)[1:(min(11, n_distinct(prob)))])
    ) %>% 
    mutate(val = prob / max(c(-Inf, prob), na.rm = T)) %>% 
    ungroup() %>% 
    arrange(bactGroup, desc(prob)) %>% ungroup() %>% 
    mutate(species = str_replace(species, "sp\\d+", "sp.")) %>% 
    left_join(bactGenomeSize %>% select(genus, species, size) %>% distinct(),
              by = c("genus", 'species')) %>% 
    mutate(runId = {{runId}},
           pipelineId = toProcess$pipelineId[sampleIndex]) %>% 
    left_join(
      myData %>% group_by(taxid, geneId) %>% 
        filter(fullPath == max(fullPath)) %>% dplyr::slice(1) %>% 
        group_by(taxid) %>% 
        summarise(
          depth = sum(KC)/sum(LN),
          pathScore = sum(pathScore), 
          .groups = "drop"
        ),
      by = "taxid"
    ) 
  
  return(list(myClusters = myClusters, bactList = bactList))
}


#Get all the ARG that are most likely found in genomes
# allBact$taxid = as.character(allBact$taxid)

genomeARG = allBact %>% 
  group_by(geneId, plasmid) %>% 
  summarise(path = max(fullPath), .groups = "drop") %>% 
  group_by(geneId) %>% 
  mutate(perc = path / max(path)) %>% 
  filter((plasmid & perc < 0.5) | !any(plasmid))

genomeARG = allBact %>% filter(geneId %in% genomeARG$geneId) 


#Check if there are any, proceed accordingly
if(nrow(genomeARG) > 0){
  
  #Adjust the bact presence across ARG
  genomeARG = adjustBact(genomeARG)
  
  #Get the clusters and bact list
  AMRclusters = genomeARG$myClusters %>%
    select(-primary) %>% filter(val > 0) %>% 
    left_join(allBact %>% 
                select(taxid, geneId, pathScore), by = c("taxid", "geneId")) %>% 
    mutate(origin = "genome")
  
  AMRclusters = AMRclusters %>% 
    left_join(
      AMRclusters %>% 
        group_by(geneId, pathScore) %>% summarise(contenders = n(), .groups = "drop") %>% 
        group_by(geneId) %>% 
        arrange(desc(pathScore)) %>% 
        mutate(contenders = cumsum(contenders)) %>% ungroup(),
      by = c("geneId", "pathScore")
    )
  
  bactList = genomeARG$bactList %>% group_by(bactGroup) %>% 
    mutate(x = depth, depth = weighted.mean(depth, prob)) %>% ungroup()
}


#Check for likely plasmids
if(n_distinct(genomeARG$myClusters$geneId) < n_distinct(allBact$geneId)){
  
  #Check if any of the remaining ARG match a genome ARG species
  checkPlasmid = allBact %>% 
    filter(!geneId %in% genomeARG$myClusters$geneId) %>% 
    left_join(
      bactList %>% 
        # filter(val == 1) %>% #maybe remove?
        select(taxid, bactGroup, prob, genomeDepth = depth,
               genomePathscore = pathScore, val),
      by = "taxid"
    ) %>% left_join(
      AMRclusters %>% 
        # filter(val == 1) %>%
        select(taxid, geneId) %>% 
        left_join(genesDetected %>% 
                    # mutate(cover = ifelse(type == "fragmentsOnly", cover2, cover)) %>% 
                    select(geneId, cover, startDepth, type),
                  by = "geneId") %>% 
        group_by(taxid) %>% 
        summarise(
          genomeCover = weighted.mean(cover, cover),
          genomeType = case_when(
            all(type == "noFragments") ~ "noFragments",
            all(type == "fragmentsOnly") ~ "fragmentsOnly",
            TRUE ~ "mixed"
          ), .groups = "drop"),
      by = "taxid"
    ) 
  
  checkPlasmid =  checkPlasmid %>% left_join(
    checkPlasmid%>% 
      group_by(geneId, pathScore) %>% summarise(contenders = n(), .groups = "drop") %>% 
      group_by(geneId) %>% 
      arrange(desc(pathScore)) %>% 
      mutate(contenders = cumsum(contenders)) %>% ungroup(),
    by = c("geneId", "pathScore")
  ) 
  
  #Pathscore should be at in the top range
  checkPlasmid = checkPlasmid %>% 
    group_by(geneId) %>% 
    # filter(pathScore >= quantile(pathScore, 0.75)) %>%
    mutate(top = pathScore == max(pathScore)) %>% 
    ungroup() %>% 
    #Give a slight edge to genome matches over plasmid
    mutate(pathScore = ifelse(!plasmid & !is.na(bactGroup), 
                              pathScore * 1.01, pathScore))
  
  
  myPlasmids = list()
  temp = checkPlasmid
  #FInd the best matching plasmids
  while(nrow(temp > 0)){
    x = temp %>%
      filter(top) %>% 
      group_by(accession, genus, species) %>% 
      mutate(x = sum(fullPath)) %>% ungroup() %>% 
      filter(x == max(x))
    
    
    myPlasmids[[length(myPlasmids)+1]] = x
    
    
    temp = temp %>% filter(!geneId %in% x$geneId) 
  }
  
  #Link genome taxid to the plasmids and asssign if match criteria
  # sapply(myPlasmids, function(x) "2518" %in%  x$geneId)
  # x = myPlasmids[[1]]
  temp = map_df(myPlasmids, function(x){
    
    test = checkPlasmid %>% filter(geneId %in% x$geneId) %>% 
      group_by(taxid, geneId) %>% filter(fullPath == max(fullPath)) %>% 
      dplyr:: slice(1) %>% 
      group_by(taxid) %>% mutate(x = sum(pathScore)) %>% 
      ungroup() %>% mutate(x = x / max(x)) %>% 
      group_by(taxid, geneId) %>% 
      summarise(pathScore = sum(pathScore), depth = mean(depth),
                across(c(bactGroup:genomeType), max), .groups = "drop") %>%
      mutate(pathScore = ifelse(is.na(bactGroup), pathScore, pathScore*1.05)) %>% 
      group_by(geneId) %>% 
      filter(depth >= 0.75 * genomeDepth | is.na(genomeDepth) |
               abs(depth - genomeDepth) < 10) %>%
      filter(pathScore == max(pathScore))
      # arrange(desc(pathScore)) %>% dplyr::slice(1)
    
    
  }, .id = "plasmidGroup") %>% 
    group_by(plasmidGroup) %>% 
    filter((bactGroup == bactGroup[!is.na(bactGroup)][1] )| all(is.na(bactGroup)))
  
  #Build data frame with the annotated plasmids
  myPlasmids = map_df(myPlasmids, function(x) x, .id = "plasmidGroup") %>% 
    select(plasmidGroup, geneId, bestAccession = accession) %>% 
    left_join(temp %>% select(plasmidGroup, taxid), by = "plasmidGroup") %>% 
    left_join(checkPlasmid, by = c("geneId", "taxid")) 
  
}

#Certainty can also be derived from complexity of gfa, the more segments and links
#the more complex and thus the more possible confusion
# score 1 = two paths, best scenaria
# score >1 = one path, actually fragment ...
# score < 1, the smaller the more paths, the more complex
pathComplexity = pathData %>% group_by(geneId) %>% 
  summarise(pathCompl = round(2/n_distinct(pathId), 3)) %>% 
  mutate(geneId = as.integer(geneId))

#Link plasmids to genomes if possible
test = bind_rows(
  allBact %>% 
    left_join(
      AMRclusters %>% filter(val == 1) %>% 
        select(taxid, geneId, bactGroup),
      by = c("taxid", "geneId")
    ) %>% filter(!is.na(bactGroup)) %>% 
    group_by(geneId, taxid) %>% filter(fullPath == max(fullPath)) %>% 
    dplyr::slice(1) %>% ungroup() %>% distinct(),
  
  #Plasmid
  allBact %>% 
    left_join(
      myPlasmids %>% distinct() %>% 
        # group_by(plasmidGroup, geneId) %>% 
        # filter(pathScore >= 0.75*max(pathScore, na.rm = T)) %>% 
        # filter(all(!bactGroup %in% genomeARG$myClusters$bactGroup) | 
        #          bactGroup %in% genomeARG$myClusters$bactGroup) %>% 
        group_by(bactGroup, geneId) %>% 
        filter(pathScore == max(pathScore)) %>% dplyr::slice(1) %>% 
        select(accession, geneId, bactGroup, plasmidGroup) %>% distinct(),
      by = c("accession", "geneId")
    ) %>% filter(!is.na(plasmidGroup)) 
  
) %>% left_join(pathComplexity, by = "geneId") %>% 
  select(bactGroup, genus, species, geneId, gene, subtype, type, pathCompl, 
         fullPath, everything()) %>% 
  arrange(bactGroup, gene, subtype)

sample

# myPlasmids %>% group_by(plasmidGroup, geneId) %>% 
#   filter(fullPath == max(fullPath)) %>% 
#   dplyr::slice(1) %>% ungroup() %>% 
#   select(geneId, taxid, accession, plasmidGroup) %>% 
#   left_join(bactList %>% select(taxid, bactGroup),
#             by = "taxid") %>% select(-taxid)