#!/bin/bash

baseFolder=$(realpath -- "$(dirname -- "$0")")

echo "1) Check dependencies..."
#Check is R is installed
test=`grep -oP "rscript\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
test=`command -v $test`
if [ -z "$test" ]; then 
	echo -e "\e[91mR does not seem to be installed.\n If it is, set the path to 'Rscript' in the settings file\n"\
	 $baseFolder/settings.txt"\e[0m"
	exit 1;
fi;

#Check if the correct R packages are installed
Rscript $baseFolder/dataAndScripts/setup.R

#Check if bbmap is installed or the reformat.sh script can be reached
test=`grep -oP "reformat\s*=\s*\K([^\s]+)" $baseFolder/settings.txt`
test=`command -v $test`
if [ -z "$test" ]; then 
	echo -e "\e[91mThe bbmap package does not seem to be installed as a system application\n"\
	"If you have unzipped the package in a custom folder,\n update the path to the 'reformat.sh' script in the settings file\n"\
	$baseFolder/settings.txt"\e[0m"
	exit 1;
fi;

echo -e "\e[32m All dependnecies seem to have been installed\e[0m\n"

#Test the whole pipeline
#-----------------------

echo -e "2) Test mixing metagenome...\n"
#Create input file
cat $baseFolder/dataAndScripts/testData/input.csv | awk '{gsub(/~/,"'$baseFolder'")}1' > \
	$baseFolder/dataAndScripts/testData/testInput.csv

#Run mixMultiple.sh
$baseFolder/mixMultiple.sh -f -i $baseFolder/dataAndScripts/testData/testInput.csv \
	-o $baseFolder/dataAndScripts/testData/testOutput.fastq.gz


