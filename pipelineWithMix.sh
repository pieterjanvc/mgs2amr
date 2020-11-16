#!/bin/bash
#BSUB -W 00:30
#BSUB -n 4
#BSUB -M 20000
#BSUB -o /data/aplab/ARG_PJ/clusterLogs/%J_pipelineWithMix.out
#BSUB -e /data/aplab/ARG_PJ/clusterLogs/%J_pipelineWithMix.err

# module load R/4.0.2
# module load sqlite3
source ~/pushover.sh

baseFolder=/mnt/meta2amrData/meta2amr

#Save error to temp file to it can be both displayed to user and put in DB
touch $baseFolder/dataAndScripts/lastError
exec 2>$baseFolder/dataAndScripts/lastError

#When error occurs, notify and exit
err_report() {

    #Use the line number where error occured and the saved error message
    errMsg=`cat $baseFolder/dataAndScripts/lastError` 
	
	#Insert into DB (make sure quoting is all right)
	errMsg=$(sed 's/'\''/&&/g' <<< "$errMsg")
    updateDBwhenError "$runId" "ERROR LINE $1: $errMsg"
	
	#Report error to stdout too 
	echo -e "\n\e[91m--- ERROR LINE $1 ---\n"
	echo -n "$errMsg"
	echo -e "\e[0m"
	
	pushMessage "Error in line $1" "pipelineWithMix"
	
	exit 1;
}
trap 'err_report ${LINENO}' ERR

updateDBwhenError() {
	#Update the DB
    $sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$2'
	WHERE runId = $1"
}

pushMessage "Start Pipeline" "pipelineWithMix"

#Folders
mixTempFolder=/mnt/meta2amrData/meta2amr/temp
meta2amrTempFolder=/mnt/meta2amrData/meta2amr/temp
genomesFolder=/mnt/meta2amrData/cchmcData/mixedMetagenomes #Folder to save the mixed files to or get existing ones
bgFolder=/mnt/meta2amrData/cchmcData/normalMetagenomes/ #Background metagenomes folder
iFolder=/mnt/meta2amrData/ncbi/sra #Isolates folder

#Variables
fileName=testMixWithBG #Set to a specific filename if starting from existing mixed file (rest below is ignored then)
bgName="H7Heidi4.fastq.gz" #Name of the background metagenome
# iNames="SRR2976831,SRR4025843" #Name(s) of the isolate(s). Comma separate (no space) if multiple
# relAb="0.07,0.05" #Sum = 1 if no bgMeta or <1 otherwise. Comma separate (no space) if multiple
iNames="SRR4025843" #Name(s) of the isolate(s). Comma separate (no space) if multiple
relAb="0.1" #Sum = 1 if no bgMeta or <1 otherwise. Comma separate (no space) if multiple

sqlite3=`grep -oP "sqlite3\s*=\s*\K(.*)" $baseFolder/settings.txt`

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (scriptName,start,status) \
	values('pipelineWithMix.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
echo -e "\n\e[32m"`date "+%T"`" - Starting pipeline with mix ...\e[0m";

if [ ! -f $genomesFolder/$fileName.fastq.gz ]; then

	# STEP 1 - Generate the mixing input file
	#-----------------------------------------
	
	printf "%s - Generate the mixing input file ... " `date "+%T"`

	#Generate the csv file that will serve as input
	Rscript $baseFolder/dataAndScripts/generateMixInput.R \
		"$genomesFolder" "$bgFolder" "$bgName" "$iFolder" \
		"$iNames" "$relAb" "$fileName"
		
	printf "done\n"
		
		
	# STEP 2 - Create the mixed file
	#-------------------------------
	
	$baseFolder/mixMultiple.sh \
		-i $genomesFolder/$fileName.csv \
		-o $genomesFolder/$fileName.fastq.gz \
		-t $mixTempFolder
	
else
	
	echo -e `date "+%T"` "-" $fileName.fastq.gz "already exists, not creating new one"
fi


# STEP 3 - meta2amr
# ------------------
$baseFolder/meta2amr.sh \
	-i $genomesFolder/$fileName.fastq.gz \
	-o $meta2amrTempFolder \
    -t $meta2amrTempFolder \
	-v 1

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
	
# cd /database/ncbi/
# myFasta=/data/aplab/ARG_PJ/aim2/mixedSample_1581975575_toBlast.fasta
# outputFolder=/data/aplab/ARG_PJ/aim2

# module load blast/2.10.0
# blastn -db nt -query $myFasta -task megablast -num_threads 1 -evalue 0.0001 -out "/data/aplab/ARG_PJ/aim2/1581975575_alignment.out" \
	# -outfmt "6 qseqid qlen sacc evalue bitscore length pident nident gaps gapopen staxids scomnames sskingdoms stitle qcovs qcovhsp qcovus" 
	
# makeblastdb -in $myFasta -parse_seqids -blastdb_version 5 -out /data/aplab/ARG_PJ/aim2/testBlastDB -title "testBlastDB" -dbtype nucl

# blastn -db /data/aplab/ARG_PJ/aim2/testBlastDB -query $myFasta -task megablast -num_threads 1 -evalue 0.0001 -out "/data/aplab/ARG_PJ/aim2/1581975575_selfAlign.out" \
	# -outfmt "6 qseqid qlen sacc evalue bitscore length pident nident gaps gapopen staxids scomnames sskingdoms stitle qcovs qcovhsp qcovus" 

pushMessage "Finished" "pipelineWithMix"
