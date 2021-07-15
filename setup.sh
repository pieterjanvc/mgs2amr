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
    $sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$2'
	WHERE runId = $1"
}

while getopts ":ht" opt; do
  case $opt in
	h) echo -e "\n"
	   awk '/--- SETUP.SH ---/,/-- END SETUP.SH ---/' $baseFolder/readme.txt
	   echo -e "\n"
	   exit
    ;;
	t) runTests=true
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
sqlite3=`grep -oP "sqlite3\s*=\s*\K(.*)" $baseFolder/settings.txt`
testTool=`command -v $sqlite3`
if [ -z "$testTool" ]; then 
    message="SQLite 3 does not seem to be installed.\n If it is, set the path to 'sqlite3' in the settings file"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	exit 1;
fi;
echo -e " - SQLite 3 is present"

#Check the meta2amr database and create if needed
if [ ! -f "$baseFolder/dataAndScripts/meta2amr.db" ]; then
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" -cmd \
	".read $baseFolder/dataAndScripts/createMeta2amrDB.sql" \
	".mode csv" ".import $baseFolder/dataAndScripts/argTable.csv ARG"
	echo -e " - No meta2amr database found, a new database was created"
else 
	echo -e " - The meta2amr database is present"
fi

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (pipelineId,scriptName,start,status) \
	values(0,'checkSetup.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
#Check if R is installed
Rscript=`grep -oP "rscript\s*=\s*\K(.*)" $baseFolder/settings.txt`
if [ -z `command -v $Rscript` ]; then 
    message="R does not seem to be installed.\n If it is, set the path to 'Rscript' in the settings file"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "R does not seem to be installed"
	exit 1;
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),1,'R installed')"

#Check if the correct R packages are installed
$Rscript $baseFolder/dataAndScripts/setup.R
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.R',$(date '+%s'),2,'R packages installed')"
echo -e " - R and dependent packages are present"

#Check if usearch is installed or the reformat.sh script can be reached
testTool=`grep -oP "usearch\s*=\s*\K(.*)" $baseFolder/settings.txt`
if [ -z `command -v $testTool` ]; then 
	echo -e "\e[91mThe usearch package does not seem to be installed as a system application\n"\
	"If you installed it in a custom folder,\n update the path to the usearch script in the settings file\n"\
	$baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "The usearch package does not seem to be installed"
	exit 1;
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),4,'usearch installed')"
echo -e " - usearch is present"

#Check if MetaCherchant.sh can be reached
testTool=`grep -oP "metacherchant\s*=\s*\K(.*)" $baseFolder/settings.txt`
if [ -z `command -v $testTool` ]; then 
	echo -e "\e[91mThe MetaCherchant script cannot be found\n"\
	"Update the path to the script in the settings file\n"\
	$baseFolder/settings.txt"\e[0m"
	updateDBwhenError "$runId" "MetaCherchant script not found"
	exit 1;
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),5,'MetaCherchant present')"
echo -e " - MetaCherchant is present"

#Check if BLASTn is either a local tool or link to a remote service
blastPath=`grep -oP "localBlastBlastn\s*=\s*\K(.*)" $baseFolder/settings.txt`
if [ -z `command -v $blastPath` ]; then 
	blastPath=`grep -oP "remoteBlastBlastn\s*=\s*\K(.*)" $baseFolder/settings.txt`
	curl -s --head $blastPath| head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
	if [ "$?" != 0 ]; then
		message="No valid path to either a local or remote BLASTn service is set\n Set the path to 'blastn' in the settings file"
		echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
		updateDBwhenError "$runId" "No valid path to either a local or remote BLASTn"
		exit 1;
	fi
	message="local BLASTn present - "
	echo -e " - Local BLASTn instance not found: localBlast.sh can NOT be used"
else
	message="local BLASTn NOT present - "
	echo -e " - Local BLASTn instance present: localBlast.sh can be used"
fi

blastPath=`grep -oP "remoteBlastBlastn\s*=\s*\K(.*)" $baseFolder/settings.txt`
curl -s --head $blastPath| head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
if [ "$?" != 0 ]; then
	message=$message"remote BLASTn present"
	echo -e " - Remote BLASTn service present: remoteBlast.sh can be used"
else
	message=$message"remote BLASTn NOT present"
	echo -e " - Remote BLASTn service present: remoteBlast.sh can be used"
fi

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),6,'$message')"

#Check if pigz is installed else use gzip (slower but same result)
if [ -z `command -v pigz` ]; then 
	message="pigz not installed. gzip used instead"
else
	message="pigz present"
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),7,'$message')"
echo -e " - $message"


#Check if pandoc version 1.12.3+ is present
if [ -z `command -v pandoc` ]; then 
	message="pandoc is NOT present. HTML reports cannot be generated"
else
	message="pandoc present"
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),7,'$message')"
echo -e " - $message"

echo -e "   ... finished\n"
finalMessage=" All dependencies seem to be present\n"

if [ "$runTests" == "true" ]; then
	#STEP 2 - Test the whole pipeline
	#---------------------------------
	printf "2) Test mixing metagenome... "

	#Create input file
	cat $baseFolder/dataAndScripts/testData/input.csv | awk '{gsub(/~/,"'$baseFolder'")}1' > \
		$baseFolder/dataAndScripts/testData/testInput.csv

	#Run mixMultiple.sh
	$baseFolder/mixMultiple.sh -f \
		-i $baseFolder/dataAndScripts/testData/testInput.csv \
		-o $baseFolder/dataAndScripts/testData/testOutput.fastq.gz \
		-v FALSE

	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
		"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
		VALUES($runId,'setup.sh',$(date '+%s'),8,'mixMultiple test succesful')"
	printf "done\n\n"
	
	finalMessage="$finalMessage  Mixing test successful\n"
fi

#Finish script
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"
	
echo -e `date "+%T"`" - Setup check finished succesfully\n \e[32m$finalMessage\e[0m"
