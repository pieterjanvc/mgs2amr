########## META2AMR ##########
##############################
   Developed by PJ Van Camp

--- SETUP.SH ---
Run the setup.sh script to verify all dependencies and to test the pipeline.

Arguments [h|t]
 -h Read the help documentation
 -t Run pipeline tests with dummy data (will take some time)

The following software needs to be installed:
- SQLite3 
  * Precompiled 32-bit version: https://www.sqlite.org/download.html
  * Precompiled 64-bit version: https://github.com/boramalper/sqlite3-x64/releases
- R version 4.0+
  * Packages: dplyr(1.0+), stringr, tidyr, jsonlite
  * Precompiled versions: https://www.r-project.org/ 
- bbmap 
  * https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
- usearch
  * https://www.drive5.com/usearch/download.html
- MetaCherchant
  * https://github.com/ctlab/metacherchant
- BLASTn: two options available
  * Use a local blastn tool 
    Make sure to download the nucleotide (nt) database or a subset containing all bacteria
   https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download    
  * Use a cloud instance of BLAST
    Example: https://blast.ncbi.nlm.nih.gov/Blast.cgi

Optional software
- pigz
  * If installed, zipping files will be faster on multi-core machines compared to gzip
   
IMPORTANT: Update the paths to all dependencies in the 'settings.txt' file 
 if they are not in the default PATH

-- END SETUP.SH ---


--- MIXMULTIPLE.SH ---
Mix multiple isolate WGS files together to create artificial metagenomes 

Arguments [h|i|o|r|m|t|f|v]
 -h Read the help documentation
 -i The input file (.csv) containing all samples to be mixed
 -o The location to save the output file. Filename should end with a fastq.gz extension
 -r (optional) The max number of reads the mixed file should contain. 
     By default, the number of reads in the background file is chosen.
     If no background file is present, the limit is the sum of 
	 the fractions needed from each isolate file.
 -m (optional) TRUE or FALSE. 
     Generate a meta-data JSON file in the same folder as the output file.
     Default can be changed in the settings.txt file
 -t (optional) The location of the temp folder (files removed once completed). 
     Default can be changed in the settings file
 -f (optional) If set, force overwriting an exisiting output file
 -v (optional) TRUE or FALSE. Progress is posted to stdout when TRUE.
     Default can be changed in the settings.txt file


Input file details --
This is a .csv file with the following columns
 - type: either I for isolate file or B for background file
   * Minimum of 2 I files if no B file
   * Max 1 B file and 1 or more I files
 - sampleName: name for the different input files (can be blank)
 - relativeAbundance: relative abundance of the file in the final metagenome.
    The sum of all must be 1.0
 - readFile: full path to the first read file (fastq.gz format)
 - readFile2: full path to the second read file (fastq.gz format)
    Leave empty in case of 1 interleaved data file
 - Any other columns will be ignored, but put in the meta-data JSON file if generated

-- END MIXMULTIPLE.SH ---


--- META2AMR.SH ---
Evaluate a metagenome for pathogenic bacterial species and their AMR

Arguments [h|i|o|v]
 -h Read the help documentation
 -i The input file in fastq or fastq.gz format
 -s (optional) The secondary input file in case of pair-end reads
 -o The folder to save the results. A sub-folder will be created
 -n The name of the output sub-folder and prefix for other filenames
 -t (optional) The location of the temp folder (files removed once completed). 
     Default = meta2amr/temp unless 'meta2amrTemp' is set in the settings file
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout
    1: General progress is posted to stdout
    2: All available details are posted to stdout    

Special case when resuming an existing pipeline
 -p Provide only a pipelineId (found in file called pipelineId in temp/output folder). 
     Use if an error occurred and you like to resume from a specific temp file
	 All other parameters but verbose will be ignored
	
-- END META2AMR.SH ---


--- LOCALBLAST.SH ---
Run local BLASTn for files in the meta2amr database awaiting alignment.

This script is separate from the main pipeline as local BLASTn searches require
a large amount of memory and thus can be run in a separate process

Arguments [b|d|h|r|v]
 -h Read the help documentation
 -b (optional) The link to the remote blastn API
	 Default = value from 'localBlastBlastn' in the settings file
 -d (optional) Set the BLAST database to use (should be nt or custom set of bacteria)
	 Default = value from 'localBlastDB' in the settings file
 -p (optional) Run the BLAST search only for these pipelineIds (e.g. "1,5,20")
     Default = all pipelineIds in the blastSubmissions table of the meta2amr database 
	 that have not been blasted yet
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout (except errors)
    1: General progress is posted to stdout
	
-- END LOCALBLAST.SH ---


--- REMOTEBLAST.SH ---
Run a remote BLASTn service for files in the meta2amr database awaiting alignment
The script will keep running until all searches have been completed or
a timeout has been reached.

Arguments [b|e|f|r|t|v]
 -h Read the help documentation
 -b (optional) The path to the local blastn module
	 Default = value from 'remoteBlastBlastn' in the settings file
 -e (optional) Entrez query to filter the nt database
	 Default = value from 'remoteBlastEntrez' in the settings file
	 Note: changing this might have unforeseen effects on the correct species calling
 -f (optional) The frequency (sec) with which the the script looks for updates on searches 
	 Default = value from 'remoteBlastCheckFreq' in the settings file
 -p (optional) Run the BLAST search only for these pipelineIds (e.g. "1,5,20")
     Default = all pipelineIds in the blastSubmissions table of the meta2amr database 
	 that have not been blasted yet
 -t (optional) The timeout (sec) after which the the script stops looking new results
	 The script will end before the timeout if all searches have been completed
	 Default = value from 'remoteBlastTimeout' in the settings file
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout (except errors)
    1: General progress is posted to stdout
	
-- END REMOTEBLAST.SH ---
