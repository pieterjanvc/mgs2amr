#####################################################
# ALTERNATIVE PIPELINE MODEL USED FOR BENCHMARKING 
####################################################

##### VARIABLES TO SUPPLY #####
# i: isolate SRA ID to use for mix-in by SEQ2MGS
# b: background metagenome SRA ID to use by SEQ2MGS
# r: relative abundance of isolate in background
# p: continue from existing pipelineId

while getopts ":i:b:r:p:" opt; do
  case $opt in
	i) isolate="${OPTARG}"
    ;;
	b) background="${OPTARG}"
    ;;
	r) ra="${OPTARG}"
    ;;
	p) pipelineId="${OPTARG}"
    ;;
    \?) echo "Unknown argument provided"
	    exit
	;;
  esac  
done


#Tools required by pipeline
#--------------------------
#spades/3.15.5
#python3
#sratoolkit/3.0.3
#R/4.2.3
#bbmap/38.96.0
#sqlite3/3.33.0
#pigz/2.6.0
#blast/2.13.0
#diamond/2.1.8
#SEQ2MGS (see previous publication)

#NOTE: make sure all tools are accessible by $PATH

#File and folder locations to set once
baseFolder="revisions"
database="metaspadesPipe.db"
export BLASTDB="/database/blast/db"

#Notify on error
set -e
function errNote {
  echo "ERROR revisionPipeline"
}
trap errNote EXIT

fmtDate () { echo $(date +'%Y-%m-%d %T' 2> /dev/null); }

#START PIPELINE
#--------------

echo "Start revision pipeline $(fmtDate date)"

#Create database if needed
if [ ! -f "$database" ]; then
	sqlite3 "$database" ".read $baseFolder/createMetaspadesPipeDB.sql"
fi

if [ "$pipelineId" == "" ]; then
	#Create new entry
	name=$isolate\_$background\_$ra
	statusCode=1
	
	pipelineId=$(sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO pipeline (statusCode,isolate,background,ra,created,modified) \
	values($statusCode,'$isolate','$background',$ra,'$(date '+%F %T')','$(date '+%F %T')'); \
	SELECT pipelineId FROM pipeline WHERE pipelineId = last_insert_rowid();")
	
	name=$name\_$pipelineId
	folder=$baseFolder/data/$name
	sqlite3 -cmd ".timeout 30000" $database \
		"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START SEQ2MGS'); 
		UPDATE pipeline SET name = '$name', folder = '$folder' WHERE pipelineId = $pipelineId"	
	
else
  
  #Continue with existing pipelineId
	pipelineId=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT pipelineId FROM pipeline WHERE pipelineId = $pipelineId;")
	
	if [ "$pipelineId" == "" ]; then 
		pushMessage "ERROR: pipelinedId not found" "CCHMC"
		exit 1
	fi
	
	name=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT name FROM pipeline WHERE pipelineId = $pipelineId;")
	folder=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT folder FROM pipeline WHERE pipelineId = $pipelineId;")
	statusCode=$(sqlite3 -cmd ".timeout 30000" $database \
	"SELECT statusCode FROM pipeline WHERE pipelineId = $pipelineId;")	
	
fi

#STEP 1 - MIX FILE - SEQ2MGS
#--------------------------
if [ "$statusCode" == 1 ]; then
	mkdir $folder
	
	echo "Run SEQ2MGS - Mix $isolate in $background with $ra RA"
	#Create seq2mgs input file
	echo "type,relativeAbundance,getFromSRA" > $folder/mixInput.csv
	echo "I,$ra,$isolate" >> $folder/mixInput.csv
	echo "B,,$background" >> $folder/mixInput.csv
	#Run seq2mgs
	seq2mgs.sh -i $folder/mixInput.csv -o $folder/$name.fastq.gz -t /tmp
	
	statusCode=2
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START metaSPAdes'); 
	UPDATE pipeline SET statusCode = 2, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip SEQ2MGS, already run"
fi

#STEP 2 - ASSEMBLE - metaspades
#------------------------------
if [ "$statusCode" == 2 ]; then
	echo "Run metaspades"
	metaspades.py --12 $folder/$name.fastq.gz -o $folder -m 64 -t 8
	
	statusCode=3
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START DIAMOND'); 
	UPDATE pipeline SET statusCode = 3, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip metaspades, already run"
fi

#STEP 3 - FIND ARG - DIAMOND
#--------------------------
if [ "$statusCode" == 3 ]; then
	echo "Run DIAMOND"
	#Only keep contigs >= 250 bp
	reformat.sh fastaminlen=250 in=$folder/scaffolds.fasta out=$folder/queryseq.fasta overwrite=true

	#Run DIAMOND
	diamond blastx -d $baseFolder/diamond/diamondDB -q $folder/queryseq.fasta -o $folder/diamondOutput.csv \
		-f 6 qseqid sseqid qlen slen pident nident length evalue bitscore qstart qend
		
	statusCode=4
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START BlastPrep'); 
	UPDATE pipeline SET statusCode = 4, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip DIAMOND, already run"
fi


#STEP 4 - BLASTN PREP - R
#--------------------------
if [ "$statusCode" == 4 ]; then
	echo "Run blast prep"
	Rscript $baseFolder/revisions_blastPrep.R "$baseFolder" "$pipelineId" "$database"
		
	statusCode=5
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START BLASTn'); 
	UPDATE pipeline SET statusCode = 5, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip blast prep, already run"
fi

#STEP 5 - ALIGN TO nt - blastn
#------------------------------
if [ "$statusCode" == 5 ]; then
	echo "Run blastn"
	blastn -db /database/blast/db/nt -query $folder/toBlast.fasta -task megablast -word_size 64 \
		-max_target_seqs 500 -max_hsps 3 -taxidlist "$baseFolder/../dataAndScripts/bact.txids" -num_threads 8 \
		-perc_identity 75 -outfmt "6 qseqid sallacc staxids sscinames salltitles qlen slen qstart qend sstart \
			send bitscore score length pident nident qcovs qcovhsp" -out $folder/blastResults.csv
		
	statusCode=6
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'START Annotation'); 
	UPDATE pipeline SET statusCode = 6, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip BLASTn, already run"
fi

#STEP 6 - Filter BLASTn output - R
#---------------------------------
if [ "$statusCode" == 6 ]; then
	echo "Run blast output filter"
	Rscript $baseFolder/revisions_filterBlast.R "$pipelineId" "$database"
		
	statusCode=7
	sqlite3 -cmd ".timeout 30000" $database \
	"INSERT INTO log (pipelineId,timeStamp,logCode,logMsg) values ($pipelineId,'$(fmtDate)',$statusCode,'Finished pipeline'); 
	UPDATE pipeline SET statusCode = 7, modified = '$(date '+%F %T')' WHERE pipelineId = $pipelineId"
else
	echo "Skip BLASTn output filtering, everything has been run already!"
fi

#PIPELINE FINISHED!
echo "Finished revision pipeline $(fmtDate date)
