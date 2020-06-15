#!/bin/bash

# source /data/aplab/ARG_PJ/pushover.sh

#When error occurs, notify and exit
err_report() {
    #pushMessage "ERROR line $1 - pipeline exit" "Pipeline AIM 2"
	echo "ERROR line $1 - pipeline exit" "Pipeline AIM 2"
	exit
}
trap 'err_report $LINENO' ERR

baseFolder=$(realpath -- "$(dirname -- "$0")/..")

#Options when script is run
while getopts ":hio" opt; do
  case $opt in
	h) cat $baseFolder/readme.txt
	   exit
    ;;	
	i) inputFile="${OPTARG%/}"
    ;;
	o) outputFile="${OPTARG%/}"
    ;;	
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done

#Check all the input arguments
if [ -z ${inputFile+x} ]; then 
	echo -e "\n\e[91mNo input file found, type mixMultiple -h for more info\e[0m"; exit 1; 
fi

if [ -z ${outputFile+x} ]; then 
	echo -e "\n\e[91mNo output file found, type mixMultiple -h for more info\e[0m"; exit 1; 
fi



#pushMessage "mixMultiple" "DONE"
