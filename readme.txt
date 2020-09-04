
###### META2AMR ######
######################
Developed by PJ Van Camp

--- SETUP.SH ---
Run the setup.sh script to verify all dependencies and to test the pipeline

The following software needs to be installed:
- R (version 4.0+ recommended)
  * Packages: dplyr, stringr, tidyr, jsonlite
- bbmap (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
  * If unzipped in custom folder, set the path to the reformat.sh script in the settings.txt
- MetaCherchant (https://github.com/ctlab/metacherchant)

Update the paths to dependencies in the 'settings.txt' file if needed

-- END SETUP.SH --


--- MIXMULTIPLE.SH ---
Mix multiple isolate WGS files together to create artificial metagenomes 

Arguments [h|i|o|r|m|t|f|v]
 -h Read the help documentation
 -i The input file (.csv) containing all samples to be mixed
 -o The location to save the output file. Filename should end with a fastq.gz extension
 -r The max number of reads the mixed file should contain. 
     By default, the number of reads in the background file is chosen.
     If no background file is present, the limit is the sum of the fractions needed from each isolate file.
 -m TRUE or FALSE. Generate a meta data JSON file in the same folder as the output file.
     Default = TRUE, but can be changed in the settings.txt file
 -t The location of the temp folder (files removed once completed). 
     Default = meta2amr/temp unless 'mixMultipleTemp' is set in the setting file
 -f If set, force overwriting an exisiting output file
 -v TRUE or FALSE. Progress is posted to stdout when TRUE.
     Default = TRUE, but can be changed in the settings.txt file
    

Input file details --
This is a .csv file with the following columns
 - type: either I for isolate file or B for background file
   * Minimum of 2 I files if no B file
   * Max 1 B file and 1 or more I files
 - sampleName: optional name for the different input files
 - relativeAbundance: relative abundance of the file in the final metagenmome.
    The sum of all must be 1.0
 - readFile: full path to the first read file (fastq.gz format)
 - readFile2: full path to the second read file (fastq.gz format)
    Leave empty in case of 1 interleaved data file
 - Any other columns will be ignored, but put in the meta data JSON file if generated

-- END MIXMULTIPLE.SH --


--- META2AMR.SH ---
Evaluate a metagenome for pathogenic bacterial species and their AMR

Arguments [h|i|o|v]
 -h Read the help documentation
 -i The input file in fastq.gz format
 -o The folder to save the results. A subfolder will be created
 -v 0: Nothing is written to stdout
    1: General progress is posted to stdout
    2: All available details are posted to stdout
    Default: meta2amrVerbose = 1, but can be changed in the settings file
	
-- END META2AMR.SH --
