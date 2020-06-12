#!/bin/bash

source /data/aplab/ARG_PJ/pushover.sh

#When error occurs, notify and exit
err_report() {
    #pushMessage "ERROR line $1 - pipeline exit" "Pipeline AIM 2"
	echo "ERROR line $1 - pipeline exit" "Pipeline AIM 2"
	exit
}
trap 'err_report $LINENO' ERR

function join { local IFS="$1"; shift; echo "$*"; } #custom function to convert Bash array to string that can be passed to R

#Options when script is run
while getopts ":f:i:s:p:b:t:l:o:" opt; do
  case $opt in
	f) isolatesFolder="${OPTARG%/}"
    ;;
    i)  set -f # disable glob
        IFS=' ' # split on space characters
        iSamples=($OPTARG)
    ;;
	s)  set -f # disable glob
        IFS=' ' # split on space characters
        iSamples2=($OPTARG)
    ;;
	p)  set -f # disable glob
        IFS=' ' # split on space characters
        procentIsolateInMix=($OPTARG)
    ;;
	b) mSample="$OPTARG"
    ;;
	t) mSample2="$OPTARG"
    ;;
	l) readLimit="$OPTARG"
    ;;
	o) outputFolder="${OPTARG%/}"
    ;;	  
    \?) echo "Unknown argument provided"
	;;
  esac  
done

#Check all the input arguments
if [ -z ${iSamples+x} ]; then 
	echo -e "\nNo isolate samples chosen. Set -i (sampleName1 sampleName2) \n"; exit 1; 
elif [ ${#iSamples[@]} -eq 1 ] && [ -z ${mSample+x} ]; then
	echo -e "\nOnly 1 sample set. Either choose at least one other sample or background \n"; exit 1; 
fi
if [ -z ${isolatesFolder+x} ]; then 
	echo -e "\nIsolates folder missing. Set -f /pathToFolder/ \n"; exit 1
else
	for file in $iSamples
	do
		if [ ! -f $isolatesFolder/$file.fastq.gz ]; then
			echo -e "\nThe following file does not exist: $isolatesFolder/$file.fastq.gz \n"; exit 1
		fi
	done
	
	if [ ! -z ${iSamples2+x} ]; then 
		for file in $iSamples2
		do
			if [ ! -f $isolatesFolder/$file.fastq.gz ]; then
				echo -e "\nThe following file does not exist: $isolatesFolder/$file.fastq.gz \n"; exit 1
			fi
		done
		
		if [ ${#iSamples[@]} -ne ${#iSamples2[@]} ]; then
		echo -e "\nThe total number of samples in -i does not equal that of -s \n"; 
		exit 1;
		fi
		
		in2='$isolatesFolder/${iSamples2[$x]}\.fastq.gz'
	else
		unset in2
	fi
fi

if [ -z ${procentIsolateInMix+x} ]; then 
	echo -e "\nNo procentIsolateInMix.  Example: set -p \"0.3 0.7\" \n"; 
	exit 1;
elif [ ${#procentIsolateInMix[@]} -ne ${#iSamples[@]} ]; then
	echo -e "\nThe number of samples does not equal the number of procentIsolateInMix \n"; 
	exit 1;
fi
if [ -z ${mSample+x} ]; then 	
	if [ $(echo "$( IFS="+"; bc <<< "${procentIsolateInMix[*]}" ) < 1.0" |bc -l) == 1 ]; then
	    echo -e "\nThe mix-in percentages do not add up to 1. This is required when there is no background file \n";
		exit 1
	elif [ $(echo "$( IFS="+"; bc <<< "${procentIsolateInMix[*]}" ) > 1.0" |bc -l) == 1 ]; then
	    echo -e "\nThe sum of mix-in percentages is greater than 1 \n";
		exit 1
	fi
	
	mSample=none
else
	if [ ! -f $mSample ]; then
			echo -e "\nThis background file does not exist: $mSample \n"; exit 1
	fi
		
	if [ $(echo "$( IFS="+"; bc <<< "${procentIsolateInMix[*]}" ) > 1.0" |bc -l) == 1 ]; then
	    echo -e "\nThe sum of mix-in percentages is greater than 1 \n";
		exit 1
	fi	
fi
if [ -z ${readLimit+x} ]; then readLimit=0; fi;
if [ -z ${outputFolder+x} ]; then outputFolder=$isolatesFolder; fi;
if [ ! -z ${mSample2+x} ]; then
	if [ ! -f $mSample2 ]; then
			echo -e "\nThe following file does not exist: $mSample2 \n"; exit 1
	fi
fi

echo $in2

#FOLDERS
baseFolder=$(realpath -- "$(dirname -- "$0")/..")

#Getting ready
nIsolates=${#iSamples[@]}
mixedSample=mixedSample_`date +%s`
mkdir -p $baseFolder/temp/$mixedSample
tempFolder=$baseFolder/temp/$mixedSample

#Write meta-data to the output folder
touch $outputFolder/$mixedSample.txt
echo -e "creationDate:" `date` \
	"\nbaseMetagenome:" $mSample \
	"\nisolates:" ${iSamples[@]} \
	"\nrelativeAbundance:" ${procentIsolateInMix[@]} >> $outputFolder/$mixedSample.txt

touch $baseFolder/log
echo -e "\n--" `date` "-- START mixing metagenomes \n" >> $baseFolder/log

#Change variable to string
procentIsolateInMix=`join , ${procentIsolateInMix[@]}`
echo `date +"%T"` mixing $iSample in $mSample with abundance $procentIsolateInMix >> $baseFolder/log

#Get the number of reads in the files
#We stored the number of reads when BBmap interleaved for metagenome (both end-reads are counted separately!) 
if [ $mSample != none ]; then 
	nReadsM=`grep -oP "(?<="$mSample", )(\d+$)" /data/aplab/ARG_PJ/data/haslamData/normalMetagenomes/readCounts.csv`
else
	nReadsM=0
fi

#We count the number of lines in one isolate read file, the number of reads in the interlveaved = nLines / 4 for singele file, 2 for pair
nReadsI=()
for x in ${iSamples[@]}; do
  nReadsI+=( $((`zcat $isolatesFolder/$x\.fastq.gz | wc -l` / `if [ -z ${iSamples2+x} ]; then echo 4; else echo 2; fi;`)) )
done
nReadsI=`join , ${nReadsI[@]}`

#Use R script to get the percentage of reads needed to get the relative abundances required
module load bbmap
module load R/4.0.0
sampleRate=(`Rscript $baseFolder/dataAndScripts/mixInAmount.R $nReadsM $nReadsI $procentIsolateInMix $readLimit`)

echo `date +"%T"` nReadsM, nReadsI, mixin, procedure: $nReadsM $nReadsI $procentIsolateInMix $procedure >> $baseFolder/log
echo `date +"%T"` sample rate: `join , ${sampleRate[@]}` >> $baseFolder/log

nReadsImixedIn=()
for x in ${!iSamples[@]}; do
     
	echo -e "\nSTART file $x"
	#Save the number of reads added
	totalReads=()
	
	#If exactly the whole file is needed, just copy paste (or interleave)
	if [ $(echo "${sampleRate[$x]} == 1.0" | bc -l) == 1 ]; then
	
		if [ -z ${iSamples2+x} ]; then
			cp $isolatesFolder/${iSamples[$x]}\.fastq.gz $tempFolder/sample$x.fastq.gz			
		else
		    reformat.sh in=$isolatesFolder/${iSamples[$x]}\.fastq.gz in2=`eval echo $in2` \
			out=$tempFolder/sample$x.fastq.gz
		fi
		echo "Copy whole file"
		totalReads+=${nReadsI[$x]}
	
	else #In case less or more than the whole file is needed
	
		#Get the integer value of the sample rate
		duplicates=`echo ${sampleRate[$x]} | grep -oP "^[0-9]+"`
		
		#If integer > 0, add a full interleaved file, but change the ID so every read stays unique when combined later
		if [ $duplicates -ne 0 ]; then		
			
			#Add as many duplicates (with different id) as needed. e.g. 1.7 = 1; 2.03 = 2
			for y in `seq $duplicates`; do
				echo "sample $x, duplicate $y"			
				reformat.sh in1=$isolatesFolder/${iSamples[$x]}\.fastq.gz \
				in2=`eval echo $in2` \
				out=stdout.fastq | awk 'NR % 4 == 1{sub(/@/,"'@$y'_",$0);print;next}NR % 2 == 1{print "+";next}{print}' | gzip -c \
				> $tempFolder/sample$x\_$y.fastq.gz
				totalReads+=${nReadsI[$x]}
			done

		fi
		
		#Add the random fraction of the sample rate: e.g. 1.7 = 0.7; 2.03 = 0.03
		fraction=`echo ${sampleRate[$x]} | sed "s/.*\./0./g"`
		
		if [ $(echo "$fraction == 0" | bc -l) == 0 ]; then	
			echo "add fraction sample $x: $fraction"
			echo $(eval echo $in2)
			totalReads+=$(reformat.sh --samplerate=$fraction in1=$isolatesFolder/${iSamples[$x]}\.fastq.gz \
			in2=$(eval echo $in2) out=$tempFolder/sample$x.fastq.gz \
			2> >( grep -oP "\d+(?= reads \()"))
		fi			
    fi	
	
	nReadsImixedIn+=`printf '%s\n' "${totalReads[@]}" | awk '{s+=$1} END {print s}'`
	nReadsImixedIn+=", "
done

#Add the random fraction of the background metagenome if set
if [ $mSample != none ]; then 
	nReadsImixedIn+=`reformat.sh --samplerate=${sampleRate[-1]} in=$mSample in2=$mSample2 \
	out=$tempFolder/background.fastq.gz 2> >( grep -oP "\d+(?= reads \()")`
fi

echo $tempFolder
#Merge all files
find $tempFolder -name '*.fastq.gz' -exec cat {} \; > $outputFolder/$mixedSample.fastq.gz
#cat $tempFolder/*.fastq.gz > $outputFolder/$mixedSample.fastq.gz

module unload bbmap

echo -e "Number of isolate reads mixed in:" $nReadsImixedIn >> $outputFolder/$mixedSample.txt
echo -e "As percentage of reads in original file:" `join , ${sampleRate[@]::$nIsolates}` >> $outputFolder/$mixedSample.txt

#Clear the temp folder
#rm -r $tempFolder

echo `date +"%T"` mix-in FINISHED >> $baseFolder/log

#pushMessage "mixMultiple" "DONE"
