#*************************
# ---- MGS2AMR Setup ----
#************************

args = commandArgs(trailingOnly = TRUE)
folder = file.path(normalizePath(args[[1]]), "dataAndScripts")

# Check assets
url = "https://github.com/pieterjanvc/mgs2amr/releases/download/assets_v1.0.0/mgs2amr_assets.tar.gz"

files = c(
  file.path(folder, "bact.txids"),
  file.path(folder, "ARG.fasta"),
  file.path(folder, "bactGenomeSize.csv"),
  file.path(folder, "testData", "testOutput.fastq.gz"),
  file.path(folder, "..", "tools", "metacherchant.sh"),
  file.path(folder, "..", "tools", "gfaTools_v0.9.0-beta.tar.gz")
)

if (!all(sapply(files, file.exists))) {
  cat(" - Downloading missinng assets ... ")
  
  dir.create(file.path(folder, "..", "tools"), showWarnings = F)
  dir.create(file.path(folder, "testData"), showWarnings = F)

  # Download the data
  download.file(url, file.path(folder, "assets.tar.gz"), mode = "wb", quiet = T)
  untar(file.path(folder, "assets.tar.gz"), exdir = folder)
  unlink(file.path(folder, "assets.tar.gz"))

  dir.create(file.path(folder, "..", "tools"), showWarnings = F)
  # Put in correct folders
  x <- file.rename(
    file.path(folder, "testOutput.fastq.gz"),
    file.path(folder, "testData", "testOutput.fastq.gz")
  )
  x <- file.rename(
    file.path(folder, "metacherchant.sh"),
    file.path(folder, "..", "tools", "metacherchant.sh")
  )
  x <- file.rename(
    file.path(folder, "gfaTools_v0.9.0-beta.tar.gz"),
    file.path(folder, "..", "tools", "gfaTools_v0.9.0-beta.tar.gz")
  )
  cat("done\n")
} else {
	cat(" - All assets present\n")
}

#Check R version
if(as.integer(R.version$major) < 4){
  warning(paste0("R version 4.0 or higher is recommended.\n",
      "Currently ", R.version$version.string, " is the default.\n",
      " If your system has multiple R versions and the default < 4.0,\n",
      " set the default 'Rscript' in $PATH to version 4+.\n",
      " In case of unexpected errors, update R to the latest version\n"))
}

#Check R packages
packages = c("tidyverse", "RSQLite", "igraph", "foreach", "doParallel", "seqinr")
installed = packages %in% installed.packages()[,1]

if(!all(installed)){
  stop(paste("The following R packages are not installed:\n",
             paste(packages[!installed], collapse = ", ")))
}

if(!stringr::str_detect(as.character(packageVersion("dplyr")), "^1")){
  stop("The dplyr package needs to be version 1.0+")
}

# Check gfaTools
if(! "gfaTools" %in% installed.packages()[,1]){
  packagePath <- normalizePath(file.path(folder, "..", "tools", "gfaTools_v0.9.0-beta.tar.gz"),mustWork = T)
  stop("The 'gfaTools' package needs to be installed manually (not on CRAN).",
       "It can be found in the mgs2amr/tools/ folder.\n", "Open R and run:",
       sprintf("\n install.packages(\"%s\")"), packagePath)
}
