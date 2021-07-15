#!/bin/bash

baseFolder=$(realpath -- "$(dirname -- "$0")")
sqlite3=`grep -oP "sqlite3\s*=\s*\K(.*)" $baseFolder/settings.txt`

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
while getopts ":hi:s:o:n:t:fv:p:m:" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- META2AMR.SH ---/,/-- END META2AMR.SH ---/' $baseFolder/readme.txt  
	   echo -e "\n"
	   exit
    ;;	
	i) inputFile1="${OPTARG}"
    ;;
	s) inputFile2="${OPTARG}"
    ;;
	o) outputFolder=`realpath "${OPTARG%/}"`
    ;;
	n) outputName="${OPTARG}"
    ;;
	t) tempFolder=`realpath "${OPTARG%/}"`
    ;;
	f) forceOverwrite=TRUE
    ;;
	v) verbose="${OPTARG}"
    ;;
	p) pipelineId="${OPTARG}"
    ;;
	m) memory="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

exec 2>$baseFolder/dataAndScripts/lastError

#Check the arguments
if [ -z ${forceOverwrite+x} ]; then 
	forceOverwrite=`grep -oP "meta2amrForceOverwrite\s*=\s*\K(.*)" $baseFolder/settings.txt`
fi

if [ -z ${memory+x} ]; then 
	# memory=$(expr $(free -h | grep -oP "Mem:\s+[^\s]+\s+[^\s]+\s+\K([\d\.]+)") - 4)
	memory="32G"
fi

if [ -z ${pipelineId+x} ]; then 	

	if [ ! -f "$inputFile1" ]; then 
		echo -e "\n\e[91minputFile1 does not exist.\n Use -i to specify one, or read the help file (-h)\e[0m"; exit 1;	
	elif [ ! -z ${inputFile2+x} ] && [ ! -f "$inputFile2" ]; then
		echo -e "\n\e[91minputFile2 does not exist.\n Use -s to specify one, or read the help file (-h)\e[0m"; exit 1; 	
	fi

	if [ -z ${outputFolder+x} ]; then 
		echo -e "\n\e[91mNo output folder specified.\n Use -o to specifiy one or read the help file (-h)\e[0m"; exit 1;
	elif [ ! -d $outputFolder ]; then	
		echo -e "\n\e[91mThe output directory does not exist\e[0m"; exit 1;
	fi

	if [ -z ${tempFolder+x} ]; then 
		tempFolder=`grep -oP "meta2amrTemp\s*=\s*\K(.*)" $baseFolder/settings.txt`
		tempFolder=${tempFolder%/}
		if [ -z ${tempFolder} ]; then
			tempFolder=$baseFolder/temp #Use default temp if none assigned
		elif [ ! -d `dirname $tempFolder` ]; then	
			echo -e "\n\e[91mThe default temp directory set in the settings file does not exist\e[0m"; exit 1;
		fi
	elif [ ! -d `dirname $tempFolder` ]; then	
		echo -e "\n\e[91mThe temp directory does not exist\e[0m"; exit 1;
	fi
	
	# if [ $forceOverwrite == FALSE ] &  [ -d `dirname $outputFolder` ]; then 
		# echo -e "\n\e[91mThe temp directory does not exist\e[0m"; exit 1;
	# fi
	
	if [ -z ${outputName+x} ]; then 
		rand=`shuf -zer -n5  {A..Z} {a..z} {0..9}`
		outputName=meta2amrPipeline
		
		while [ -d $outputFolder/$outputName$rand ]; do
			rand=`shuf -zer -n5  {A..Z} {a..z} {0..9}`
		done
		
		outputName=$outputName\_$rand	
    elif [ ! $(grep -E "^\\w+$" <<< $outputName) ]; then
		echo -e "\n\e[91mThe (-n) outputName can only consist of alphanumeric and underscore characters\e[0m"; exit 1;
	elif [ -d $outputFolder/$outputName ]; then
		echo -e "\n\e[91mA folder with name $outputName already exists in the output folder." \
		"Use -n to specifiy a unique name or read the help file (-h)\e[0m"; exit 1;
	fi
	
	tempName=$outputName\_`date '+%s'`
	
else

	#Get the first runId
    firstRunId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT s.runId 
	FROM scriptArguments as s, scriptUse as p \
	WHERE p.pipelineId = $pipelineId AND p.runId = s.runId AND \
	      s.scriptName = 'meta2amr.sh' AND argument = 'inputFile1'")
	
	if [ -z "$firstRunId" ]; then
		firstRunId=0
	fi
	
	tempFolder=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'meta2amr.sh' AND argument = 'tempFolder' AND runId = $firstRunId")
	
	tempName=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'meta2amr.sh' AND argument = 'tempName' AND runId = $firstRunId")
	
	if [ ! -f "$tempFolder/$tempName/pipelineId" ]; then
		echo -e "\n\e[91mThe the pipelineId provided does not point to a valid folder \e[0m"; exit 1; 
	fi
	
    #In case of a previous runId, load all the arguments from the database
    inputFile1=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'meta2amr.sh' AND argument = 'inputFile1' AND runId = $firstRunId")
	
	inputFile2=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'meta2amr.sh' AND argument = 'inputFile2' AND runId = $firstRunId")
	
	outputFolder=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'meta2amr.sh' AND argument = 'outputFolder' AND runId = $firstRunId")

	MCsuccess=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"SELECT logId FROM logs
	WHERE runId IN (SELECT runId FROM scriptUse WHERE pipelineId == $pipelineId) AND \
	actionId = 2")
	
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "meta2amrVerbose\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif ! grep -qE "^(0|1|2)$" <<< $verbose ; then
echo verbose = ${verbose[@]}
	echo -e "\n\e[91mThe verbose option (-v) needs to be 0, 1 or 2\n Read the help file (-h) for more info\e[0m"; exit 1; 
fi

if [ $verbose != "0" ]; then echo -e "\n\e[32m"`date "+%T"`" - Inputs correct, starting pipeline ...\e[0m"; fi;


#Register the start of the script in the DB
if [ -z ${pipelineId+x} ]; then
	#Generate the next pipelineId
	pipelineId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO pipeline (name,tempFolder,outputFolder,statusCode,statusMessage,startTimestamp,modifiedTimestamp) \
	values('$outputName','$tempFolder/$tempName','$outputFolder/$outputName',1,'Pipeline started','$(date '+%F %T')','$(date '+%F %T')'); \
	SELECT pipelineId FROM pipeline WHERE pipelineId = last_insert_rowid()")
	
    runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values($pipelineId,'meta2amr.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
	scriptArgs="($runId,'meta2amr.sh','inputFile1', '$inputFile1'),
	($runId,'meta2amr.sh','inputFile2', '$inputFile2'),
	($runId,'meta2amr.sh','outputFolder', '$outputFolder'),
	($runId,'meta2amr.sh','outputName', '$outputName'),
	($runId,'meta2amr.sh','tempFolder', '$tempFolder'),
	($runId,'meta2amr.sh','tempName', '$tempName'),"	
else
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE pipeline \
	SET modifiedTimestamp = '$(date '+%F %T')' \
	WHERE pipelineId == $pipelineId"
	
	runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values($pipelineId, 'meta2amr.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
	pointerTopipelineId="($runId,'meta2amr.sh','pipelineId', '$pipelineId'),"
	
fi

#Save the arguments with which the script was run
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES $pointerTopipelineId $scriptArgs
	($runId,'meta2amr.sh','forceOverwrite', '$forceOverwrite'),
	($runId,'meta2amr.sh','verbose', '$verbose')"


#--- PART 1 METACHERCHANT (MC) ---
#---------------------------------
	
#Only run if MetaCherchant has not been run before for this sample (i.e. the temp folder has MC data)
echo "*****************************"
echo "--- STEP 1: MetaCherchant ---"
echo "*****************************"
if [ -z "$MCsuccess" ]; then 
	
	if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Start MetaCherchant ..."; fi;
	
	#Generate temp folders
  mkdir -p $tempFolder/$tempName
	mkdir -p $tempFolder/$tempName/metacherchant_logs
	rm -rf $tempFolder/$tempName/metacherchant_logs/*
	echo $pipelineId > $tempFolder/$tempName/pipelineId
    
	#MC generates a lot of output, we ignore this but when verbose = 2
	if [ $verbose != 2 ] ; then
		exec 2>/dev/null #Ignore error stream (info) of metacherchant
	else
		exec 2>&1 #Rederect stream to stdout
    fi
	
	metacherchant=`grep -oP "metacherchant\s*=\s*\K(.*)" $baseFolder/settings.txt`
	
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),1,'Start MetaCherchant')"
	
	inputFile="$inputFile1 $inputFile2"
	$metacherchant --tool environment-finder \
		--k 31 \
		--coverage=5 \
		--reads $inputFile \
		--seq $baseFolder/dataAndScripts/ARG.fasta \
		--maxradius 2000 \
		--output $tempFolder/$tempName \
		--work-dir $tempFolder/$tempName/metacherchant_logs \
		--maxkmers=100000 \
		--bothdirs=False \
		--chunklength=250 \
		-m $memory
		
    $sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),2,'Finished MetaCherchant')"	
	
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE pipeline 
	SET statusCode = 2, statusMessage = 'Finished MetaCherchant', modifiedTimestamp = '$(date '+%F %T')'
	WHERE pipelineId == $pipelineId"
		
	exec 2>$baseFolder/dataAndScripts/lastError #Restore output redirection
	
	if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Finished MetaCherchant"; fi;
	

else

	if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Skip MetaCherchant, already done"; fi;
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),3,'Skip MetaCherchant, already done')"
	
	touch $tempFolder/$tempName/pipelineId
	
fi


#--- PART 2 BLAST PREP ---
#-------------------------
echo -e "\n"
echo "***********************************"
echo "--- STEP 2: BLASTn Preparations ---"
echo "***********************************"

if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Start BLAST preparations  ..."; fi;

#Get paths from the settings
Rscript=`grep -oP "rscript\s*=\s*\K(.*)" $baseFolder/settings.txt`

#Set any option to modify the script
scriptOptions=(keepAllMetacherchantData maxPathDist minBlastLength trimLength clusterIdentidy forceRedo, maxStep)
scriptValues=(FALSE 1500 250 250 0.99 "$forceOverwrite" 0)

if [ "${scriptValues[6]}" != 0 ]; then 
	echo "  WARNING: the blast prep is limited to step ${scriptValues[6]}" 
fi

# declare -A scriptOptions
# scriptOptions[keepAllMetacherchantData]=FALSE #keepAllMetacherchantData, if FALSE, only the GFA data is kept after merging (folders are removed)
# scriptOptions[maxPathDist]=5000 #maxPathDist: Distance from ARG to crop the GFA file (reduces blast search)
# scriptOptions[minBlastLength]=250 #minBlastLength: Min segment length to submit to blast
# scriptOptions[trimLength]=100 #trimLength: Loose segments smaller than this will be cut from thr GFA
# scriptOptions[clusterIdentidy]=0.95 #clusterIdentidy: The cluster identity percent used in usearch
# scriptOptions[forceRedo]=FALSE #forceRedo: if TRUE, all code is run again, even is intermediate files are available

for i in "${!scriptOptions[@]}"; 
do 
	# values=$values",('$runId','$i','${scriptOptions[$i]}')" 
	values=$values",('$runId','${scriptOptions[$i]}','${scriptValues[$i]}')"
done
values=`echo $values | sed -e 's/^,//g'`

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO blastPrepOptions (runId,option,value) VALUES $values"

$Rscript $baseFolder/dataAndScripts/blastPrep.R \
	"$baseFolder" "$tempFolder"	"$tempName" "$verbose" "$runId" "$pipelineId" \
	${scriptValues[@]}

if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Finished BLAST preparations"; fi;

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
