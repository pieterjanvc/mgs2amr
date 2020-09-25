
###### META2AMR ######
######################
Developed by PJ Van Camp

--- SETUP.SH ---
Run the setup.sh script to verify all dependencies and to test the pipeline

The following software needs to be installed:
- SQLite3 (https://www.sqlite.org/download.html)
- R version 4.0+ (https://www.r-project.org/) 
  * Packages: dplyr(1.0+), stringr, tidyr, jsonlite
- bbmap (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
  * If unzipped in custom folder, set the path to the reformat.sh script in the settings.txt
- MetaCherchant (https://github.com/ctlab/metacherchant)
- BLASTn: two options available
  * Use a local blastn tool (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    Make sure to download the nucleotide (nt) database or a subset containing all bacteria
  * Use a cloud instance of BLAST (e.g. https://blast.ncbi.nlm.nih.gov/Blast.cgi)	
   
IMPORTANT: Update the paths to all dependencies in the 'settings.txt' file if needed

-- END SETUP.SH ---


--- MIXMULTIPLE.SH ---
Mix multiple isolate WGS files together to create artificial metagenomes 

Arguments [h|i|o|r|m|t|f|v]
 -h Read the help documentation
 -i The input file (.csv) containing all samples to be mixed
 -o The location to save the output file. Filename should end with a fastq.gz extension
 -r (optional) The max number of reads the mixed file should contain. 
     By default, the number of reads in the background file is chosen.
     If no background file is present, the limit is the sum of the fractions needed from each isolate file.
 -m (optional) TRUE or FALSE. Generate a meta-data JSON file in the same folder as the output file.
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
 - relativeAbundance: relative abundance of the file in the final metagenmome.
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
 -i The input file in fastq.gz format
 -o The folder to save the results. A subfolder will be created
 -t (optional) The location of the temp folder (files removed once completed). 
     Default = meta2amr/temp unless 'meta2amrTemp' is set in the setting file
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout
    1: General progress is posted to stdout
    2: All available details are posted to stdout    

Special case when resuming previous run
 -r Provide only a previous runId (found in file called runId in temp folder). 
     Only set if an error occured and you like to resume from a specific temp file
	 All other parameters but verbose will be ignored
	
-- END META2AMR.SH ---


--- LOCALBLAST.SH ---
Run local BLASTn for files in the meta2amr database awaiting alignment

Arguments [b|d|h|r|v]
 -h Read the help documentation
 -b (optional) The path to the local blastn module
	 Default = value from 'localBlastBlastn' in the setting file
 -d (optional) Set the BLAST database to use (should be nt or custom set of bacteria)
	 Default = value from 'localBlastDB' in the setting file
 -r (optional) Run the BLAST search only for these runIds 
     Default = all runIds in the blastSubmissions table of the meta2amr database 
	 that have not been blasted yet
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout (except errors)
    1: General progress is posted to stdout
	
-- END LOCALBLAST.SH ---
