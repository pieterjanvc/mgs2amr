#!/bin/bash

#When error occurs, notify and exit
err_report() {
    errMsg=`cat $baseFolder/dataAndScripts/lastError`
    echo -e "$errMsg" #report error to stdout too
    updateDBwhenError "Error Line $1 $errMsg"
	exit
}
trap 'err_report ${LINENO}' ERR

updateDBwhenError() {
	#Update the DB
    $sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'error',
	info = '$1'
	WHERE runId = $runId"	
}

baseFolder=$(realpath -- "$(dirname -- "$0")")

while getopts ":h" opt; do
  case $opt in
	h) echo -e "\n"
	   grep -zo "\-\-\- SETUP\.SH.*\-\- END SETUP\.SH \-\-" $baseFolder/readme.txt 
	   echo -e "\n"
	   exit
    ;;	
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done


echo "1) Check dependencies..."
#Check if sqlite3 is installed
sqlite3=`grep -oP "sqlite3\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
testTool=`command -v $sqlite3`
if [ -z "$testTool" ]; then 
    message="SQLite 3 does not seem to be installed.\n If it is, set the path to 'sqlite3' in the settings file"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	exit 1;
fi;

#Check the database and create if needed
if [ ! -f "$baseFolder/dataAndScripts/meta2amr.db" ]; then
	$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	".read $baseFolder/dataAndScripts/createMeta2amrDB.sql"
fi

exec 2>$baseFolder/dataAndScripts/lastError

#Register the start of the script in the DB
runId=$($sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO scriptUse (scriptName,start,status) \
	values('checkSetup.sh','$(date '+%F %T')','running'); \
	SELECT runId FROM scriptUse WHERE runId = last_insert_rowid()")
	
#Check if R is installed
Rscript=`grep -oP "rscript\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
testTool=`command -v $Rscript`
if [ -z "$testTool" ]; then 
    message="R does not seem to be installed.\n If it is, set the path to 'Rscript' in the settings file"
	echo -e "\e[91m$message\n" $baseFolder/settings.txt"\e[0m"
	updateDBwhenError "R does not seem to be installed"
	exit 1;
fi;

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),1,'R installed')"

#Check if the correct R packages are installed
$Rscript $baseFolder/dataAndScripts/setup.R
# exec 3>&1
# errMsg=$($Rscript $baseFolder/dataAndScripts/setup.R 2>&1 1>&3)
# exec 3>&-
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.R',$(date '+%s'),2,'R packages installed')"

#Check if bbmap is installed or the reformat.sh script can be reached
test=`grep -oP "reformat\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
test=`command -v $test`
if [ -z "$test" ]; then 
	echo -e "\e[91mThe bbmap package does not seem to be installed as a system application\n"\
	"If you have unzipped the package in a custom folder,\n update the path to the 'reformat.sh' script in the settings file\n"\
	$baseFolder/settings.txt"\e[0m"
	updateDBwhenError "The bbmap package does not seem to be installed"
	exit 1;
fi;
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),3,'bbmap installed')"

echo -e "\e[32m   All dependencies seem to have been installed\e[0m\n"


#Test the whole pipeline
#-----------------------

echo -e "2) Test mixing metagenome..."
#Create input file
cat $baseFolder/dataAndScripts/testData/input.csv | awk '{gsub(/~/,"'$baseFolder'")}1' > \
	$baseFolder/dataAndScripts/testData/testInput.csv

#Run mixMultiple.sh
$baseFolder/mixMultiple.sh -f -i $baseFolder/dataAndScripts/testData/testInput.csv \
	-o $baseFolder/dataAndScripts/testData/testOutput.fastq.gz

$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"INSERT INTO logs (runId,tool,timeStamp,actionId,actionName)
	VALUES($runId,'setup.sh',$(date '+%s'),4,'mixMultiple test succesful')"

#Update the DB
$sqlite3 "$baseFolder/dataAndScripts/meta2amr.db" \
	"UPDATE scriptUse
	SET end = '$(date '+%F %T')', status = 'finished'
	WHERE runId = $runId"