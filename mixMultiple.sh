#!/bin/bash

baseFolder=$(realpath -- "$(dirname -- "$0")")
sqlite3=`grep -oP "sqlite3\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`

#Save error to temp file to it can be both displayed to user and put in DB
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
while getopts ":hi:o:t:r:m:fv:" opt; do
  case $opt in
	h) echo -e "\n"
	   grep -zo "\-\-\- MIXMULTIPLE\.SH.*\-\- END MIXMULTIPLE\.SH \-\-" $baseFolder/readme.txt 
	   echo -e "\n"
	   exit
    ;;	
	i) inputFile=`realpath "${OPTARG}"`
    ;;
	o) outputFile=`realpath "${OPTARG}"`
    ;;
	t) tempFolder=`realpath "${OPTARG}"`
    ;;
	r) readLimit="${OPTARG}"
    ;;
	m) metaData="${OPTARG}"
    ;;
	f) forceOverwrite=T
    ;;
	v) verbose="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

exec 2>$baseFolder/dataAndScripts/lastError

#Check all the input arguments
if [ -z ${inputFile+x} ]; then 
	echo -e "\n\e[91mNo input file found, type mixMultiple -h for more info\e[0m"; exit 1; 
elif [ ! -f $inputFile ]; then 
	echo -e "\n\e[91mThe specified input file was not found\e[0m"; exit 1; 
fi

if [ -z ${outputFile+x} ]; then 
	echo -e "\n\e[91mNo output file specified.\n Use -o to specifiy one or type mixMultiple -h for more info\e[0m"; exit 1;
elif [ ! -d `dirname $outputFile` ]; then	
	echo -e "\n\e[91mThe directory for the output file does not exist\e[0m"; exit 1;
elif [ -f $outputFile ] && [ -z ${forceOverwrite+x} ]; then	
	echo -e "\n\e[91mThe output file already exists.\n Use -f option to force overwrite\e[0m"; exit 1;
fi

if [ -z ${readLimit+x} ]; then 
	readLimit=0
elif [[ ! $test =~ ^[0-9]+$ ]]; then 
	echo -e "\n\e[91mThe read limit must be a positive integer\e[0m"; exit 1; 
fi

if [ -z ${tempFolder+x} ]; then 
	tempFolder=`grep -oP "mixMultipleTemp\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
	tempFolder=${tempFolder%/}
	if [ -z ${tempFolder} ]; then
		tempFolder=$baseFolder/temp #Use default temp if none assigned
	elif [ ! -d `dirname $tempFolder` ]; then	
		echo -e "\n\e[91mThe default temp directory set in the settings file does not exist\e[0m"; exit 1;
	fi
elif [ ! -d `dirname $tempFolder` ]; then	
	echo -e "\n\e[91mThe temp directory does not exist\e[0m"; exit 1;
fi

if [ -z ${metaData+x} ]; then 
	metaData=`grep -oP "mixMultipleMetaData\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
elif ! grep -qE "^(true|T|TRUE|false|F|FALSE)$" <<< $metaData ; then	
	echo -e "\n\e[91mThe metaData option (-m) needs to be either TRUE or FALSE\e[0m"; exit 1; 
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "mixMultipleVerbose\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
elif ! grep -qE "^(true|T|TRUE|false|F|FALSE)$" <<< $verbose ; then	
	echo -e "\n\e[91mThe verbose option (-v) needs to be either TRUE or FALSE\e[0m"; exit 1; 
fi

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (scriptName,start,status) \
	values('mixMultiple.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
#Save the arguments with which the script was run
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES($runId,'mixMultiple.sh','inputFile', '$inputFile'),
	($runId,'mixMultiple.sh','outputFile', '$outputFile'),
	($runId,'mixMultiple.sh','tempFolder', '$tempFolder'),
	($runId,'mixMultiple.sh','readLimit', '$readLimit'),
	($runId,'mixMultiple.sh','metaData', '$metaData'),
	($runId,'mixMultiple.sh','forceOverwrite', '$forceOverwrite'),
	($runId,'mixMultiple.sh','verbose', '$verbose')"

echo -e "\n\e[32m"`date "+%T"`" - Start mixing reads into" `basename $outputFile` "...\e[0m"

#Run the R script
rPath=`grep -oP "rscript\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
$rPath $baseFolder/dataAndScripts/mixMultiple.R \
	$baseFolder $inputFile $outputFile $readLimit \
	$metaData $verbose $tempFolder $runId

echo -e "\e[32m"`date "+%T"`" - Finished mixing reads\n\e[0m"

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
