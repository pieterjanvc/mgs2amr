
###### META2AMR ######
######################
Developed by PJ Van Camp

--- SETUP.SH ---
Run the setup.sh script to verify all dependencies

The following software needs to be installed:
- R (version 4.0+ recommended)
- bbmap (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
- MetaCherchant (https://github.com/ctlab/metacherchant)

--- MIXMULTIPLE.SH ---
Mix multiple isolate WGS files together to create artificial metagenomes 

Arguments [-i|o]
 -i The input file (.csv) containing all samples to be mixed
 -o The location to save the output file. Filename should end with a fastq.gz extension

Input file
This is a .csv file with the following columns
- type: either I for isolate file or B for background file
	* Minimum of 2 I files if no B file
	* Max 1 B file and 1 or more I files
- sampleName: optional name for the different input files
- relativeAbundance: relative abundance of the file in the final metagenmome. The sum of all must be 1.0
- readFile: full path to the read file (fastq.gz format)
- readFile2: optional. Full path to the second read file (fastq.gz format) in case of non-interleaved pair-end data
