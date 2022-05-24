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
    sqlite3 "$database" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$2'
	WHERE runId = $1"
}

while getopts ":ht:d:" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- SETUP.SH ---/,/-- END SETUP.SH ---/' $baseFolder/readme.txt
	   echo -e "\n"
	   exit
    ;;
	t) runTests=true
	;;
	d) database="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

echo -e `date "+%T"`"\e[32m - Start the setup check...\e[0m\n"

# STEP 1 - Check dependencies
#----------------------------
echo "1) Check dependencies..."

#Check if settings file exists
if [ ! -f "$baseFolder/settings.txt" ]; then
	cp $baseFolder/dataAndScripts/backupSettings.txt $baseFolder/settings.txt
fi

#Check if sqlite3 is installed
if [ -z `command -v sqlite3` ]; then 
    message="SQLite 3 does not seem to be installed.\n If it is, add 'sqlite3' folder to $PATH"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	exit 1;
fi;
echo -e " - SQLite 3 is present"

#Check the mgs2amr database and create if needed
if [ -z ${database+x} ]; then 
	database="$baseFolder/dataAndScripts/mgs2amr.db"
fi

if [ ! -f "$database" ]; then
	sqlite3 "$database" ".read $baseFolder/dataAndScripts/createMgs2amrDB.sql"
	echo -e " - No existing mgs2amr database found, a new database was created"
else 
	echo -e " - The mgs2amr database is present"
fi

#Register the start of the script in the DB
runId=$(sqlite3 "$database" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values(0,'checkSetup.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
#Check if R is installed
if [ -z `command -v Rscript` ]; then 
    message="R does not seem to be installed.\n If it is, add the 'Rscript' folder to $PATH"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "R does not seem to be installed"
	exit 1;
fi;
sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),1,'R installed')"

#Check if the correct R packages are installed
Rscript $baseFolder/dataAndScripts/setup.R
sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.R',$(date '+%s'),2,'R packages installed')"
echo -e " - R and dependent packages are present"

#Check if usearch is installed 
if [ -z `command -v usearch` ]; then 
	echo -e "\e[91mThe 'usearch' command cannot be found\n"\
	"Set 'alias usearch=/pathToFolder/usearch<version>' to allow MGS2AMR to find the tool\n"\
	$baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "The usearch package does not seem to be installed"
	exit 1;
fi;
sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),4,'usearch installed')"
echo -e " - usearch is present"

#Check if MetaCherchant.sh can be reached
if [ -z `command -v metacherchant.sh` ]; then 
	echo -e "\e[91mThe MetaCherchant script cannot be found\n"\
	"Add the 'metacherchant.sh' folder to $PATH or create and alias\n"\
	$baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "MetaCherchant script not found"
	exit 1;
fi;
sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),5,'MetaCherchant present')"
echo -e " - MetaCherchant is present"

#Check if BLASTn is present
if [ -z `command -v blastn` ]; then 
	message="local BLASTn NOT present, "
	echo -e " - BLASTn not found\n"\
	"   Make sure blastn is installed and setup correctly (run setup.sh -h for details)"
else
	message="local BLASTn present, "
	echo -e " - Local BLASTn instance present"
fi

#Check if the $BLASTDB variable is set
if [ -z "$BLASTDB" ]; then
	BLASTDB=`grep -oP "localBlastDB\s*=\s*\K(.*)" $baseFolder/settings.txt`
fi

if blastdbcmd -list "$BLASTDB" | grep -q Nucleotide ; then 
	message=$message"nt DB detected"
	echo -e " - Nucleotide database found"
else 
	message=$message"nt DB not detected"
	echo -e " - The \$BLASTDB variable is not found\n"\
	"   Set by: export BLASTDB=/path/to/ntDBfolder"
fi

sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),6,'$message')"

#Check if pigz is installed else use gzip (slower but same result)
if [ -z `command -v pigz` ]; then 
	message="pigz not found. gzip used instead"
else
	message="pigz present"
fi;
sqlite3 "$database" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),7,'$message')"
echo -e " - $message"


echo -e "   ... finished\n"
finalMessage=" All dependencies seem to be present\n"

if [ "$runTests" == "true" ]; then
	#STEP 2 - Test the whole pipeline
	#---------------------------------
	printf "2) Test the pipeline ... "

	#... input data

	#Run mgs2amr.sh
	# $baseFolder/mgs2amr.sh -f ...

	sqlite3 "$database" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'setup.sh',$(date '+%s'),8,'pipeline test succesful')"
	printf "done\n\n"
	
	finalMessage="$finalMessage  Pipeline test successful\n"
fi

#Finish script
sqlite3 "$database" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
	
echo -e `date "+%T"`" - Setup check finished succesfully\n \e[32m$finalMessage\e[0m"
