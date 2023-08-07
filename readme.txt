 __  __  ____ ____ ____     _    __  __ ____  
|  \/  |/ ___/ ___|___ \   / \  |  \/  |  _ \ 
| |\/| | |  _\___ \ __) | / _ \ | |\/| | |_) |
| |  | | |_| |___) / __/ / ___ \| |  | |  _ < 
|_|  |_|\____|____|_____/_/   \_|_|  |_|_| \_\
                                              
     MGS2AMR - Developed by PJ Van Camp
	 
--- SETUP.SH ---
Run the setup.sh script to verify all dependencies and to test the pipeline.

IMPORTANT: In order to run this pipeline, you must download and extract 
 the zip file from the release containing all necessary files and data.
 The tracked files on GitHub only contain scripts and not other data.
 
 This pipeline will have to be run in a Linux environment 

Arguments [h|t]
 -h Read the help documentation
 -t Run pipeline tests with dummy data (will take some time)

The following software needs to be installed and available in $PATH:
- SQLite3 
  * Precompiled 32-bit version: https://www.sqlite.org/download.html
  * Precompiled 64-bit version: https://github.com/boramalper/sqlite3-x64/releases
- R version 4+ (Rscript is called)
  * Packages: tidyverse, RSQLite, igraph, foreach, doParallel, seqinr
	In addition the gfaTools package, created for this pipline needs to be installed
	The gfaTools_v0.9.0-beta.tar.gz file can be found in the root folder
  * Precompiled versions: https://www.r-project.org/ 
- usearch (v11+)
  * https://www.drive5.com/usearch/download.html
  * ! Make sure to rename the file 'usearch' if version info is attached
    mgs2amr will look for 'usearch' in $PATH
- MetaCherchant
  * This tool uses an older version of MetaChechant wich is provided in the release
  * Source https://github.com/ctlab/metacherchant
- BLASTn (Part of the BLAST+ package): 
  * https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  * Make sure to download and build the nucleotide (nt) database  
  * Set the nt database folder to the BLASTDB variable and export 
	 e.g.: `export BLASTDB=/path/to/ntDBfolder`, alternatively you can set this
	 path in the settings.txt file
  * Download the taxid database and extact in the BLASTDB folder
    or use: `perl update_blastdb.pl --passive --timeout 300 --force --verbose taxdb`
	whilst in the BLASTDB folder

Optional software
- pigz
  * If installed, zipping files will be faster on multi-core machines compared to gzip
   
IMPORTANT: Make sure all dependencies are in the $PATH variable
(sqlite3, Rscript, usearch, metacherchant.sh, blastn)

Command: 
export PATH=/pathToDependencyFolder:$PATH  
	e.g.: export PATH=/opt/ncbi-blast-2.13.0+/bin:$PATH

-- END SETUP.SH ---


--- mgs2amr.sh ---
Evaluate a metagenome for pathogenic bacterial species and their AMR
The ARG tested obtained from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047
 Last update of this list: 2022-05-31

Arguments [c|d|f|h|i|j|m|n|o|p|s|v|z]

REQUIRED
 -i The input file in fastq or fastq.gz format
 -j (optional) The secondary input file in case of pair-end reads
 -o The folder to save the results. A sub-folder will be created
    Use -n to set the name of the folder and other output files

Special case when resuming an existing pipeline
 -p Provide only a pipelineId (found in file called pipelineId in output folder). 
     The tool will resume from the last completed step if all data is available.
	 All other parameters but step (-s) and verbose (-v) will be ignored.
	 Note that the pipelineId must be found in the database (see -d)

OPTIONAL
 -c Default = 4 (or lower if not available); Number of processors to use
 -d Default = mgs2amr database in the dataAndScripts folder.
    Supply any file ending in .db to use an alternative database.
	Make sure to run the setup.sh script fist to initialise new databases
 -f Force redo and ovewrite of previous runs
 -h Read the help documentation
 -m Default = 32G; Memory available (important for MetaCherchant)
    Input files of > 2Gb easily need 32+Gb of RAM. 
	If the pipeline fails at the first step, consider increasing the memory 
 -n The name of the output sub-folder and prefix for other filenames
    If not set, a random name will be generated
 -s Set the step to which the pipeline must be run
    * 1: Step 1 - MetaCherchant
	* 2: previous + Step 2 - BlastPreparations
	* 3: previous + Step 3 - Blastn
	* 4: previous + Step 4 - Annotation
 -v Default can be changed in the settings file
    0: Nothing is written to stdout
    1: General progress is posted to stdout
    2: All available details are posted to stdout
 -z Generate output as tar.gz file in addition to storing results in the database
	Default can be changed in the settings file
	TRUE or FALSE

-- END mgs2amr.sh ---


--- LOCALBLAST.SH ---
Run local BLASTn for files in the mgs2amr database awaiting alignment.

This script can be run separately from the main pipeline as BLASTn searches require
a large amount of memory (~150 GB or more). If no specific pipelineId provided,
the scrip will seach for all runs that finished step 2 (MetaCherchant and blast prep)
and run BLASTn for all of them at once.

The first BLASTn run can take a long time as the database need to be loaded into
memory, after that, the runtime is significantly recuded for other runs in the
same session.

Arguments [d|h|r|v]
 -d Default = mgs2amr database in the dataAndScripts folder.
    Supply any file ending in .db to use an alternative database.
	Make sure to run the setup.sh script fist to initialise new databases
 -h Read the help documentation
 -p (optional) Run the BLAST search only for these pipelineIds (e.g. "1,5,20")
     Default = all pipelineIds in the blastSubmissions table of the mgs2amr database 
	 that have not been blasted yet
 -v (optional) Default can be changed in the settings file
    0: Nothing is written to stdout (except errors)
    1: General progress is posted to stdout
	
-- END LOCALBLAST.SH ---

--- EXPLORER ---
The MGS2AMR explorer Shiny app can be run in R to explore the output data.

All results are saved in the mgs2amr database and can be accessed by connecting
to it through the app. Alternatively, if a zip file with results was generated using
the -z option, the output can also be loaded that way.  

-- END EXPLORER ---