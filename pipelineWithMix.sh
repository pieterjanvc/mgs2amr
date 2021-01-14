#!/bin/bash

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

#Options when script is run
while getopts ":i:b:r:f:v:" opt; do
  case $opt in
	i) iNames="${OPTARG}"
    ;;
	b) bgName="${OPTARG}"
    ;;
	r) relAb="${OPTARG}"
    ;;
	f) fileName="${OPTARG}"
    ;;
	v) verbose="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done


if [ -z ${verbose+x} ]; then 
	verbose=1
elif ! grep -qE "^(0|1|2)$" <<< $verbose; then	
	echo -e "\n\e[91mThe verbose option (-v) needs to be 0, 1, or 2\e[0m"; exit 1;
fi

#Folders
mixTempFolder=/mnt/meta2amrData/meta2amr/temp
meta2amrTempFolder=/mnt/meta2amrData/meta2amr/temp
genomesFolder=/mnt/meta2amrData/cchmcData/mixedMetagenomes #Folder to save the mixed files to or get existing ones
bgFolder=/mnt/meta2amrData/cchmcData/normalMetagenomes/ #Background metagenomes folder
iFolder=/mnt/meta2amrData/ncbi/sra #Isolates folder

# #Variables
# fileName=testMixWithBG #Set to a specific filename if starting from existing mixed file (rest below is ignored then)
# bgName="H7Heidi4.fastq.gz" #Name of the background metagenome
# # iNames="SRR2976831,SRR4025843" #Name(s) of the isolate(s). Comma separate (no space) if multiple
# # relAb="0.07,0.05" #Sum = 1 if no bgMeta or <1 otherwise. Comma separate (no space) if multiple
# iNames="SRR4025843" #Name(s) of the isolate(s). Comma separate (no space) if multiple
# relAb="0.1" #Sum = 1 if no bgMeta or <1 otherwise. Comma separate (no space) if multiple

sqlite3=`grep -oP "sqlite3\s*=\s*\K(.*)" $baseFolder/settings.txt`

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (scriptName,start,status) \
	values('pipelineWithMix.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
#Save the arguments with which the script was run
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES($runId,'pipelineWithMix.sh','iNames', '$iNames'),
	($runId,'pipelineWithMix.sh','bgName', '$bgName'),
	($runId,'pipelineWithMix.sh','relAb', '$relAb'),
	($runId,'pipelineWithMix.sh','fileName', '$fileName'),
	($runId,'pipelineWithMix.sh','verbose', '$verbose')"
	
	
echo -e "\n\e[32m"`date "+%T"`" - Starting pipeline with mix ...\e[0m";

# STEP 1 - Download data if needed
#---------------------------------
iNameList=($(grep -oP [^,]+ <<< $test))
for mySRR in "${iNameList[@]}"
do
	if -z $(ls $iFolder | grep "$mySRR"); then
		#Download data from NCBI SRA
		echo -e "\n"`date "+%T"`" - Start downloading seq data for $mySRR ... "
		/opt/sratoolkit.2.10.8/bin/fasterq-dump $mySRR \
			-O  $outputFolder \
			-t /mnt/meta2amrData/ncbi/sra/temp
		echo -e `date "+%T"`" - downloading completed"
		
		echo -e `date "+%T"`" - Start zipping fastq files"
		#Zip the results into gz format
		find $outputFolder/$mySRR* -execdir pigz '{}' ';'
		echo -e `date "+%T"`" - finished zipping"
		
		$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),1,'Downloading $mySRR')"
		
	else
		echo -e `date "+%T"`" - $mySRR already downloaded"
		$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),2,'Skip downloading $mySRR, already present')"
	fi
done


if [ ! -f $genomesFolder/$fileName.fastq.gz ]; then

	# STEP 2 - Generate the mixing input file
	#-----------------------------------------
	
	printf "%s - Generate the mixing input file ... " `date "+%T"`

	#Generate the csv file that will serve as input
	Rscript $baseFolder/dataAndScripts/generateMixInput.R \
		"$genomesFolder" "$bgFolder" "$bgName" "$iFolder" \
		"$iNames" "$relAb" "$fileName"
		
	printf "done\n"
	
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),3,'Generated mixing input file')"
		
		
	# STEP 3 - Create the mixed file
	#-------------------------------
	
	$baseFolder/mixMultiple.sh \
		-i $genomesFolder/$fileName.csv \
		-o $genomesFolder/$fileName.fastq.gz \
		-t $mixTempFolder
		
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),4,'Mixed file created')"
	
else
	
	echo -e `date "+%T"` "-" $fileName.fastq.gz "already exists, not creating new one"
fi


# STEP 4 - meta2amr (until step 2)
# --------------------------------
$baseFolder/meta2amr.sh \
	-i $genomesFolder/$fileName.fastq.gz \
	-o $meta2amrTempFolder \
    -t $meta2amrTempFolder \
	-v $verbose

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),5,'meta2amr.sh completed until step 2')"

# STEP 5 - Delete the mix file
# ----------------------------
#We only need to keep the processed results
rm $genomesFolder/$fileName.fastq.gz

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'pipelineWithMix.sh',$(date '+%s'),6,'Mixed file deleted')"

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
