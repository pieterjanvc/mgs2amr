#!/bin/bash

#When error occurs, notify and exit
err_report() {
	echo -e "\n\e[91mERROR line $1 - meta2amr.sh\e[0m"
	exit
}
trap 'err_report $LINENO' ERR

baseFolder=$(realpath -- "$(dirname -- "$0")")

#Options when script is run
while getopts ":hi:o:fv:" opt; do
  case $opt in
	h) echo -e "\n"
	   grep -zo "\-\-\- META2AMR\.SH.*\-\- END META2AMR\.SH \-\-" $baseFolder/readme.txt 
	   echo -e "\n"
	   exit
    ;;	
	i) inputFile=`realpath "${OPTARG}"`
    ;;
	o) outputFolder=`realpath "${OPTARG%/}"`
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

#Check the arguments
if [ -z ${inputFile+x} ]; then 
	echo -e "\n\e[91mNo input file set.\n Use -i to specify one, or read the help file (-h)\e[0m"; exit 1; 
elif [ ! -f $inputFile ]; then 
	echo -e "\n\e[91mThe specified input file was not found\e[0m"; exit 1; 
elif ! grep -q "\.fastq\.gz$" <<< "mixedMetagenomes/EC_KP.fastq.gz"; then 
	echo -e "\n\e[91mThe specified input file is not in fastq.gz format\e[0m"; exit 1; 
fi

if [ -z ${outputFolder+x} ]; then 
	echo -e "\n\e[91mNo output folder specified.\n Use -o to specifiy one or read the help file (-h)\e[0m"; exit 1;
elif [ ! -d $outputFolder ]; then	
	echo -e "\n\e[91mThe output directory does not exist\e[0m"; exit 1;
fi

if [ -z ${verbose+x} ]; then 
	verbose=`grep -oP "meta2amrVerbose\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
elif ! grep -qE "^(0|1|2)$" <<< $verbose ; then
	echo -e "\n\e[91mThe verbose option (-v) needs to be 0, 1 or 2\n Read the help file (-h) for more info\e[0m"; exit 1; 
fi

if [ $verbose != "0" ]; then echo -e "\n\e[32m"`date "+%T"`" - Inputs correct, starting pipeline ...\e[0m"; fi;

#Generate temp folders
tempName=`grep -oP "\K([^\/]+)(?=.fastq.gz)" <<< $inputFile`_`date '+%s'`
mkdir $baseFolder/temp/$tempName
mkdir $baseFolder/temp/$tempName/working

if [ $verbose != "2" ]; then exec 2>/dev/null; fi;

#--- PART 1 METACHERCHANT ---
#----------------------------
if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Start MetaCherchant ..."; fi;

metacherchant=`grep -oP "metacherchant\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
$metacherchant --tool environment-finder \
	--k 31 \
	--coverage=20 \
	--reads $inputFile \
	--seq $baseFolder/dataAndScripts/ARG_06Jan2020.fasta \
	--output $baseFolder/temp/$tempName \
	--work-dir $baseFolder/temp/$tempName/working \
	--maxkmers=100000 \
	--bothdirs=False \
	--chunklength=250 \
	-m 6G

#Remove MC output files that are not needed
find $baseFolder/temp/$tempName -type d -name tsvs -delete
find $baseFolder/temp/$tempName -type f -name '*.fasta' -delete
find $baseFolder/temp/$tempName -type f -name '*.txt' -delete

if [ $verbose != "0" ]; then echo -e `date "+%T"`" -  finished"; fi;


#--- PART 2 ARG DETECTION ---
#----------------------------
if [ $verbose != "0" ]; then echo -e `date "+%T"`" - Start ARG detection ..."; fi;
