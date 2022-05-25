#!/bin/bash

baseFolder=$(realpath -- "$(dirname -- "$0")")

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
  sqlite3 $database \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$2'
	WHERE runId = $1"
}

#Options when script is run
while getopts ":hg:p:v:d:c:" opt; do
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
    d) database="${OPTARG}"
    ;;
	c) cpu="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

exec 2>$baseFolder/dataAndScripts/lastError

#Check all the input arguments
if [ -z ${database+x} ]; then 
	database=database
elif [ ! -f $database ]; then	
		echo -e "\n\e[91mThe database provided does not exist\e[0m"; exit;
fi

if [ -z ${generateReport+x} ]; then 
	generateReport=`grep -oP "annotationHTMLreport\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif [ ! $(grep -iE "^(true|false|t|f)$" <<< $generateReport) ]; then	
	echo -e "\n\e[91mThe generateHTMLReport option (-g) needs to be either true or false\e[0m"; exit 1;
fi

if [ $(grep -iE "^(true|t)$" <<< $generateReport) ]; then
	#Check if pandoc is present
	if [ -z `command -v pandoc` ]; then 
		echo -e "\n\e[91mPandoc is NOT present.\n Set -g (generateHTMLReport) to false or install Pandoc\e[0m"; 
		exit 1;
	fi;
fi

if [ -z ${cpu+x} ]; then 
	default=4
	available=`nproc --all`
	cpu=$(( available > default ? default : available ))
elif [ ! $(grep -E "^[0-9]+$" <<< $cpu) ] ; then
	echo -e "\n\e[91mThe number of processors to use needs to be a postive integer \n" \ 
	"Read the help file (-h) for more info\e[0m"; exit 1;
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "annotationVerbose\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif ! grep -qE "^(0|1|-1)$" <<< $verbose; then	
	echo -e "\n\e[91mThe verbose option (-v) needs to be either 0 or 1\e[0m"; exit 1;
fi

#Register the start of the script in the DB
runId=$(sqlite3 $database \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values(0,'annotation.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
if [ ! -z ${pipelineId+x} ]; then 
	pipelineArg="($runId,'annotation.sh','pipelineId', '$pipelineId'),"
fi
	
#Save the arguments with which the script was run
sqlite3 $database \
	"INSERT INTO scriptArguments (runId,scriptName,argument,value)
	VALUES $pipelineArg
	($runId,'annotation.sh','generateReport', '$generateReport'),	
	($runId,'annotation.sh','verbose', '$verbose')"
	
SECONDS=0

#--- PART 3 ANNOTATION ---
#-------------------------
if [ "$verbose" -gt 0 ]; then 
	echo -e "\n"
	echo "*****************************************"
	echo "--- META2AMR: ANNOTATION & PREDICTION ---"
	echo "*****************************************"
fi


#Run annotation script
Rscript $baseFolder/dataAndScripts/ARG_annotation.R \
	"$baseFolder" "$database" "$runId" "$verbose" "$pipelineId" "$generateReport" "$cpu"

#Update the DB
sqlite3 $database \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"

duration=$SECONDS

if [ "$verbose" -gt 0 ]; then 
	echo -e "\n\e[32m--- Annotation & Prediction finished successfully ---\n"\
								  "             Time elapsed: $(($duration / 60)) minutes \e[0m"
fi				
