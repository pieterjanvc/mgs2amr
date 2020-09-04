#!/bin/bash
#BSUB -W 00:10
#BSUB -n 1
#BSUB -M 8000
#BSUB -o /data/aplab/ARG_PJ/clusterLogs/%J_pipelineWithMix.out
#BSUB -e /data/aplab/ARG_PJ/clusterLogs/%J_pipelineWithMix.err

source /data/aplab/ARG_PJ/pushover.sh

#When error occurs, notify and exit
err_report() {
    pushMessage "ERROR line $1 - exit" "pipelineWithMix"
	exit
}
trap 'err_report $LINENO' ERR

# pushMessage "Start Pipeline" "pipelineWithMix"

module load R/4.0.0

#General settings
baseFolder=/data/aplab/ARG_PJ/aim2/meta2amr
tempFolder=/scratch/van9wf/pipelineTemp

mkdir -p $tempFolder/mixedMetagenomes

# STEP 1 - Generate the mixing input file
bgName="D5Heidi011" #Name of the background metagenome
iNames="ERR304796" #Name(s) of the isolate(s). Comma separate (no space) if multiple
relAb="0.2" #Sum = 1 if no bgMeta or <1 otherwise. Comma separate (no space) if multiple
fileName=testFile #Will get a timestamp extension in case name is used multiple times 

#Folders where to find the files
bgFolder=/data/aplab/ARG_PJ/data/haslamData/normalMetagenomes/
iFolder=/data/aplab/ARG_PJ/data/publicData

fileName=$fileName\_`date '+%s'`

# # STEP 1 - Generate the mixing input file
# bgName=background #Name of the background metagenome
# iNames="isolate1,isolate2" #Name(s) of the isolate(s)
# relAb="0.1,0.1" #Relative abundance of each isolate (sum must be 1 if no bgMeta or <1 otherwise)
# fileName=testFile #Will get a timestamp extension in case name is used multiple times 
# fileName=$fileName\_`date '+%s'`

# #Folders where to find the files
# bgFolder=$baseFolder/dataAndScripts/testData
# iFolder=$baseFolder/dataAndScripts/testData

#Generate the csv file that will serve as input
Rscript $baseFolder/dataAndScripts/generateMixInput.R \
	$tempFolder/mixedMetagenomes $bgFolder $bgName $iFolder \
	$iNames $relAb $fileName

# STEP 2 - Create the mixed file
$baseFolder/mixMultiple.sh \
	-i $tempFolder/mixedMetagenomes/$fileName.csv \
	-o $tempFolder/mixedMetagenomes/$fileName.fastq.gz \
	-t $tempFolder
	

# cd /database/ncbi/
# myFasta=/data/aplab/ARG_PJ/aim2/mixedSample_1581975575_toBlast.fasta
# outputFolder=/data/aplab/ARG_PJ/aim2

# module load blast/2.10.0
# blastn -db nt -query $myFasta -task megablast -num_threads 1 -evalue 0.0001 -out "/data/aplab/ARG_PJ/aim2/1581975575_alignment.out" \
	# -outfmt "6 qseqid qlen sacc evalue bitscore length pident nident gaps gapopen staxids scomnames sskingdoms stitle qcovs qcovhsp qcovus" 
	
# makeblastdb -in $myFasta -parse_seqids -blastdb_version 5 -out /data/aplab/ARG_PJ/aim2/testBlastDB -title "testBlastDB" -dbtype nucl

# blastn -db /data/aplab/ARG_PJ/aim2/testBlastDB -query $myFasta -task megablast -num_threads 1 -evalue 0.0001 -out "/data/aplab/ARG_PJ/aim2/1581975575_selfAlign.out" \
	# -outfmt "6 qseqid qlen sacc evalue bitscore length pident nident gaps gapopen staxids scomnames sskingdoms stitle qcovs qcovhsp qcovus" 

# metacherchant=/usr/local/metacherchant/1.0.0/out/metacherchant.sh
# inputFile=/data/aplab/ARG_PJ/data/test/CJ_EC.fastq.gz
# baseFolder=/data/aplab/ARG_PJ/aim2/meta2amr
# tempFolder=/scratch/van9wf/meta2amrTemp
# tempName=test


# module load R/4.0.0

# /data/aplab/ARG_PJ/aim2/meta2amr/mixMultiple.sh \
	# -i /data/aplab/ARG_PJ/data/test/SE_KP2.csv \
	# -o /data/aplab/ARG_PJ/data/test/SE_KP2.fastq.gz \
	# -f



# /data/aplab/ARG_PJ/aim2/meta2amr/meta2amr.sh \
	# -i /data/aplab/ARG_PJ/data/test/SE_KP2.fastq.gz \
	# -o /data/aplab/ARG_PJ/aim2/meta2amr/temp

pushMessage "Finished" "pipelineWithMix"
