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
	echo -e "\n\e[91m--- ERROR LINE $1 (pipelineId $pipelineId, runId $runId)---\n"
	echo -n "$errMsg"
	echo -e "\e[0m"
	
	exit 1;
}
trap 'err_report ${LINENO}' ERR

updateDBwhenError() {
	#Update the DB
    sqlite3 -cmd ".timeout 30000" $database \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$2'
	WHERE runId = $1;"
}

#Options when script is run
while getopts ":hi:j:o:n:fs:v:p:m:c:d:z:" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- mgs2amr.sh ---/,/-- END mgs2amr.sh ---/' $baseFolder/readme.txt  
	   echo -e "\n"
	   exit
    ;;	
	i) inputFile1="${OPTARG}"
    ;;
	j) inputFile2="${OPTARG}"
    ;;
	o) outputFolder=`realpath "${OPTARG%/}"`
    ;;
	n) outputName="${OPTARG}"
    ;;
	f) forceOverwrite=TRUE
    ;;
	s) step="${OPTARG}"
    ;;
	v) verbose="${OPTARG}"
    ;;
	p) pipelineId="${OPTARG}"
    ;;
	m) memory="${OPTARG}"
    ;;
	c) cpu="${OPTARG}"
    ;;
	d) database="${OPTARG}"
    ;;
	z) generateReport="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

exec 2>$baseFolder/dataAndScripts/lastError

#Check the arguments
if [ -z ${database+x} ]; then 
	database="$baseFolder/dataAndScripts/mgs2amr.db"
elif [ ! -f $database ]; then	
		echo -e "\n\e[91mThe database provided does not exist\e[0m"; exit;
fi

if [ -z ${forceOverwrite+x} ]; then 
	forceOverwrite=`grep -oP "mgs2amrForceOverwrite\s*=\s*\K(.*)" $baseFolder/settings.txt`
fi

if [ -z ${memory+x} ]; then 
	# memory=$(expr $(free -h | grep -oP "Mem:\s+[^\s]+\s+[^\s]+\s+\K([\d\.]+)") - 4)
	memory="32G"
fi

if [ -z ${cpu+x} ]; then 
	default=4
	available=`nproc --all`
	cpu=$(( available > default ? default : available ))
elif [ ! $(grep -E "^[0-9]+$" <<< $cpu) ] ; then
	echo -e "\n\e[91mThe number of processors to use needs to be a postive integer \n" \ 
	"Read the help file (-h) for more info\e[0m"; exit 1;
fi

if [ -z ${step+x} ]; then 
	step=`grep -oP "mgs2amrStep\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif [ ! $(grep -E "^(1|2|3|4)$" <<< $step) ] ; then
	echo -e "\n\e[91mThe step option (-s) needs to be between 1-4\n Read the help file (-h) for more info\e[0m"; exit 1;
fi

if [ ! $(grep -E "^(1|2|3|4)$" <<< $step) ]; then 
	echo -e "\n\e[91mThe step value (-s) must be between 1 - 4\e[0m"; exit 1;
fi

#Check if usearch is installed 
if [ -z `command -v usearch` ]; then 
	echo -e "\e[91mThe 'usearch' command cannot be found\n"\
	"Check the documentation and run setup.sh to verify\e[0m"
	exit 1;
else
	usearch=`which usearch`
fi;

if [ -z ${generateReport+x} ]; then 
	generateReport=`grep -oP "generateZipResults\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif [ ! $(grep -iE "^(true|false|t|f)$" <<< $generateReport) ]; then	
	echo -e "\n\e[91mThe generateZipResults option (-z) needs to be either true or false\e[0m"; exit 1;
fi

#Check if to create new pipelineId or check provided one
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
	
	if [ -z ${outputName+x} ]; then 
		rand=`shuf -zer -n5  {A..Z} {a..z} {0..9}`
		outputName=mgs2amr_
		
		while [ -d $outputFolder/$outputName$rand ]; do
			rand=`shuf -zer -n5  {A..Z} {a..z} {0..9}`
		done
		
		outputName=$outputName$rand	
    elif [ ! $(grep -iP "^[a-z][\w+\-\.]*\w$" <<< $outputName) ]; then
		echo -e "\n\e[91mThe (-n) outputName must be as follows:\n"\
				" - start with letter\n"\
				" - zero or more alphanumeric, '.', '_' or '-' characters\n"\
		        " - end with alphanumeric character\e[0m"; exit 1;
	elif [ -d $outputFolder/$outputName ]; then
		echo -e "\n\e[91mThe chosen output folder $outputFolder/$outputName already exists." \
		"Set -n to a unique name or read the help file (-h)\e[0m"; exit 1;
	fi
	
else

	#Get the first runId
    firstRunId=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT s.runId 
	FROM scriptArguments as s, scriptUse as p \
	WHERE p.pipelineId = $pipelineId AND p.runId = s.runId AND \
	      s.scriptName = 'mgs2amr.sh' AND argument = 'inputFile1';")
	
	if [ -z "$firstRunId" ]; then
		firstRunId=0
	fi
	
	outputFolder=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'mgs2amr.sh' AND argument = 'outputFolder' AND runId = $firstRunId;")
	
	outputName=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'mgs2amr.sh' AND argument = 'outputName' AND runId = $firstRunId;")
	
	if [ ! -f "$outputFolder/$outputName/pipelineId" ]; then
	  echo "$outputFolder/$outputName"
		echo -e "\n\e[91mThe the pipelineId provided does not point to a valid folder \e[0m"; exit 1; 
	fi
	
    #In case of a previous runId, load all the arguments from the database
  inputFile1=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'mgs2amr.sh' AND argument = 'inputFile1' AND runId = $firstRunId;")
	
	inputFile2=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'mgs2amr.sh' AND argument = 'inputFile2' AND runId = $firstRunId;")
	
	outputFolder=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT value FROM scriptArguments \
	WHERE scriptName = 'mgs2amr.sh' AND argument = 'outputFolder' AND runId = $firstRunId;")

	MCsuccess=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT logId FROM logs
	WHERE runId IN (SELECT runId FROM scriptUse WHERE pipelineId == $pipelineId) AND \
	actionId = 2;")
	
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "mgs2amrVerbose\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif [ ! $(grep -E "^(0|1|2)$" <<< $verbose) ] ; then
	echo -e "\n\e[91mThe verbose option (-v) needs to be 0, 1 or 2\n Read the help file (-h) for more info\e[0m"; exit 1; 
fi

#Set verbose to negative to let other scripts know the main echo is coming from here
verbose=-$verbose

#Register the start of the script in the DB
if [ -z ${pipelineId+x} ]; then
	#Generate the next pipelineId
	pipelineId=$(sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO pipeline (name,outputFolder,statusCode,statusMessage,startTimestamp,modifiedTimestamp) \
	values('$outputName','$outputFolder/$outputName',1,'Pipeline started','$(date '+%F %T')','$(date '+%F %T')'); \
	SELECT pipelineId FROM pipeline WHERE pipelineId = last_insert_rowid();")
	
  runId=$(sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values($pipelineId,'mgs2amr.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid();")
	
	scriptArgs="($runId,'mgs2amr.sh','inputFile1', '$inputFile1'),
	($runId,'mgs2amr.sh','inputFile2', '$inputFile2'),
	($runId,'mgs2amr.sh','outputFolder', '$outputFolder'),
	($runId,'mgs2amr.sh','outputName', '$outputName'),"	
else
	sqlite3 -cmd ".timeout 30000" $database \
	"UPDATE pipeline \
	SET modifiedTimestamp = '$(date '+%F %T')' \
	WHERE pipelineId == $pipelineId;"
	
	runId=$(sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values($pipelineId, 'mgs2amr.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid();")
	
	pointerTopipelineId="($runId,'mgs2amr.sh','pipelineId', '$pipelineId'),"
	
fi

#Save the arguments with which the script was run
sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES $pointerTopipelineId $scriptArgs
	($runId,'mgs2amr.sh','forceOverwrite', '$forceOverwrite'),
	($runId,'mgs2amr.sh','verbose', '$verbose');"

SECONDS=0

if [ "$verbose" -ne 0 ]; then
	echo -e "\n\e[32m##################################"
	echo -e "\e[32m--- MGS2AMR PIPELINE - ID $pipelineId ---"
	echo -e "\e[32m##################################\e[0m\n"
fi

if [ $step -lt 4 ]; then
	echo -e "NOTE: The pipeline is limited to step $step\n"
fi

#--- PART 1 METACHERCHANT (MC) ---
#---------------------------------
	
#Only run if MetaCherchant has not been run before for this sample (i.e. the output folder has MC data)
if [ "$verbose" -ne 0 ]; then
	echo "*****************************"
	echo "--- STEP 1: MetaCherchant ---"
	echo "*****************************"
fi

if [ -z "$MCsuccess" ]; then 
	
	if [ "$verbose" -ne 0 ]; then echo -e `date "+%T"`" - Start MetaCherchant ..."; fi;
	
	#Generate output folders
    mkdir -p $outputFolder/$outputName
	mkdir -p $outputFolder/$outputName/metacherchant_logs
	rm -rf $outputFolder/$outputName/metacherchant_logs/*
	echo $pipelineId > $outputFolder/$outputName/pipelineId
    
	#MC generates a lot of output, we ignore this but when verbose = 2
	if [ $verbose != 2 ] ; then
		exec 2>/dev/null #Ignore error stream (info) of metacherchant
	else
		exec 2>&1 #Redirect stream to stdout
    fi
	
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),1,'Start MetaCherchant');"
	
	inputFile="$inputFile1 $inputFile2"
	$baseFolder/tools/metacherchant.sh --tool environment-finder \
		--k 31 \
		--coverage=5 \
		--reads $inputFile \
		--seq $baseFolder/dataAndScripts/ARG.fasta \
		--maxradius 2000 \
		--output $outputFolder/$outputName \
		--work-dir $outputFolder/$outputName/metacherchant_logs \
		--maxkmers=100000 \
		--bothdirs=False \
		--chunklength=250 \
		-m $memory
		
    sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),2,'Finished MetaCherchant');"	
	
	sqlite3 -cmd ".timeout 30000" $database \
	"UPDATE pipeline 
	SET statusCode = 2, statusMessage = 'Finished MetaCherchant', modifiedTimestamp = '$(date '+%F %T')'
	WHERE pipelineId == $pipelineId;"
		
	exec 2>$baseFolder/dataAndScripts/lastError #Restore output redirection
	
	if [ "$verbose" -ne 0 ]; then echo -e `date "+%T"`" - Finished MetaCherchant"; fi;
	

else

	if [ "$verbose" -ne 0 ]; then echo -e `date "+%T"`" - Skip MetaCherchant, already done"; fi;
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'metacherchant.sh',$(date '+%s'),3,'Skip MetaCherchant, already done');"
	
	touch $outputFolder/$outputName/pipelineId
	
fi


#--- PART 2 BLAST PREP ---
#-------------------------
if [ $step -gt 1 ]; then
	
	if [ "$verbose" -ne 0 ]; then 
		echo -e "\n"
		echo "***********************************"
		echo "--- STEP 2: BLASTn Preparations ---"
		echo "***********************************"
	fi
	
	# if [ "$verbose" -ne 0 ]; then echo -e `date "+%T"`" - Start BLAST preparations  ..."; fi;

	#Set any option to modify the script
	scriptOptions=(keepAllMetacherchantData maxPathDist minBlastLength trimLength clusterIdentidy forceRedo, maxStep)
	scriptValues=(FALSE 2000 250 250 0.99 "$forceOverwrite" 0)

	if [ "${scriptValues[6]}" != 0 ]; then 
		echo "  WARNING: the blast prep is limited to step ${scriptValues[6]}" 
	fi

	for i in "${!scriptOptions[@]}"; 
	do 
		# values=$values",('$runId','$i','${scriptOptions[$i]}')" 
		values=$values",('$runId','${scriptOptions[$i]}','${scriptValues[$i]}')"
	done
	values=`echo $values | sed -e 's/^,//g'`

	sqlite3 -cmd ".timeout 30000" $database \
		"INSERT INTO blastPrepOptions (runId,option,value) VALUES $values"

	Rscript $baseFolder/dataAndScripts/blastPrep.R \
		"$baseFolder" "$database" "$outputFolder" "$outputName" "$verbose" "$runId" "$pipelineId" "$cpu" "$usearch" \
		${scriptValues[@]}

	if [ "$verbose" -ne 0 ]; then echo -e `date "+%T"`" - Finished BLAST preparations"; fi;
fi

if [ $step -gt 2 ]; then

	if [ "$verbose" -ne 0 ]; then 
		echo -e "\n"
		echo "**********************"
		echo "--- STEP 3: BLASTn ---"
		echo "**********************"
	fi

	if [ `command -v "blastn"` ]; then 
		$baseFolder/localBlast.sh -v "$verbose" -p "$pipelineId" -d "$database"
	else
		echo -e "\e[91mThe blastn command was not found, please make sure it's in $PATH\e[0m"
		step=2
		exit 1;
	fi
	
fi


if [ $step -gt 3 ]; then

	if [ "$verbose" -ne 0 ]; then 
		echo -e "\n"
		echo "**************************"
		echo "--- STEP 4: ANNOTATION ---"
		echo "**************************"
	fi

	$baseFolder/annotation.sh -v "$verbose" -p "$pipelineId" -d "$database" -c "$cpu" -z "$generateReport"
fi


#Update the DB
sqlite3 -cmd ".timeout 30000" $database \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
	
duration=$SECONDS

if [ $step != 4 ]; then step="at step $step"; else step=""; fi
if [ "$verbose" -ne 0 ]; then 
	echo -e "\n\e[32m--- The MGS2AMR pipeline finished successfully $step ---\n"\
						"                Time elapsed: $(($duration / 60)) minutes \e[0m"
fi
