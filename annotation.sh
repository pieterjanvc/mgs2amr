#!/bin/bash

baseFolder=$(realpath -- "$(dirname -- "$0")")
sqlite3=`grep -oP "sqlite3\s*=\s*\K(.*)" $baseFolder/settings.txt`
Rscript=`grep -oP "rscript\s*=\s*\K(.*)" $baseFolder/settings.txt`

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
while getopts ":hg:p:v:" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- ANNOTATION.SH ---/,/-- END ANNOTATION.SH ---/' $baseFolder/readme.txt
	   echo -e "\n"
	   exit
    ;;	
	g) generateReport="${OPTARG}"
	;;
	p) pipelineId="${OPTARG}"
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
if [ -z ${generateReport+x} ]; then 
	generateReport=`grep -oP "annotationHTMLreport\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif [ ! $(grep -iE "^(true|false|t|f)$" <<< $generateReport) ]; then	
	echo -e "\n\e[91mThe generateReport option (-g) needs to be either true or false\e[0m"; exit 1;
fi

if [ $(grep -iE "^(true|t)$" <<< $generateReport) ]; then
	#Check if pandoc is present
	if [ -z `command -v pandoc` ]; then 
		echo -e "\n\e[91mPandoc is NOT present. set -g to false or install Pandoc\e[0m"; 
		exit 1;
	fi;
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "annotationVerbose\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif ! grep -qE "^(0|1)$" <<< $verbose; then	
	echo -e "\n\e[91mThe verbose option (-v) needs to be either 0 or 1\e[0m"; exit 1;
fi

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values(0,'annotation.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
if [ ! -z ${pipelineId+x} ]; then 
	pipelineArg="($runId,'annotation.sh','pipelineId', '$pipelineId'),"
fi
	
#Save the arguments with which the script was run
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES $pipelineArg
	($runId,'annotation.sh','generateReport', '$generateReport'),	
	($runId,'annotation.sh','verbose', '$verbose')"

#--- PART 3 ANNOTATION ---
#-------------------------
echo -e "\n"
echo "***************************************"
echo "--- STEP 3: ANNOTATION & PREDICTION ---"
echo "***************************************"

#Run annotation script
$Rscript $baseFolder/dataAndScripts/ARG_annotation.R \
	"$baseFolder" "$runId" "$verbose" "$pipelineId" "$generateReport"

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
