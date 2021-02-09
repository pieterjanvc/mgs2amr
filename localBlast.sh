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
while getopts ":hb:d:r:v:" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- LOCALBLAST.SH ---/,/-- END LOCALBLAST.SH ---/' $baseFolder/readme.txt
	   echo -e "\n"
	   exit
    ;;	
	b) blastn=`realpath "${OPTARG}"`
    ;;
	d) blastDB=`realpath "${OPTARG}"`
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
if [ -z ${blastn+x} ]; then 
	blastn=`grep -oP "localBlastBlastn\s*=\s*\K(.*)" $baseFolder/settings.txt`
fi

if [ -z ${blastDB+x} ]; then 
	blastDB=`grep -oP "localBlastDB\s*=\s*\K(.*)" $baseFolder/settings.txt`
fi

# if [ -z ${pipelineId+x} ]; then 
	# pipelineId=0
# fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "localBlastVerbose\s*=\s*\K(.*)" $baseFolder/settings.txt`
elif ! grep -qE "^(0|1)$" <<< $verbose; then	
	echo -e "\n\e[91mThe verbose option (-v) needs to be either 0 or 1\e[0m"; exit 1;
fi

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values(0,'localBlast.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")

#Run BLASTn for all in the queue (unless runId specified)
$Rscript $baseFolder/dataAndScripts/localBlast.R \
	"$baseFolder" "$runId" "$verbose" "$blastn" "$blastDB" "$pipelineId"


#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"