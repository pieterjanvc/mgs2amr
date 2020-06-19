#!/bin/bash

#When error occurs, notify and exit
err_report() {
	echo -e "\n\e[91mERROR line $1 - mixMultiple.sh\e[0m"
	exit
}
trap 'err_report $LINENO' ERR

baseFolder=$(realpath -- "$(dirname -- "$0")")

#Options when script is run
while getopts ":hi:o:r:m:fv:" opt; do
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

echo -e "\e[32mStart mixing ...\e[0m"
#Run the R script
rPath=`grep -oP "rscript\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
$rPath $baseFolder/dataAndScripts/mixMultiple.R $baseFolder $inputFile $outputFile $readLimit $metaData $verbose
echo -e "\e[32mFinished mixing reads\n\e[0m"
