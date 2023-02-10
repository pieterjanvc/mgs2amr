#*************************
# ---- MGS2AMR Setup ----
#************************

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

if(! "gfaTools" %in% installed.packages()[,1]){
  stop("The 'gfaTools' package needs to be installed manually (not on CRAN).",
       "It can be found in the mgs2amr root folder.\n", "Open R and run:",
       "\n install.packages(\"gfaTools_<version>.tar.gz\")")
}
