#!/bin/bash
#usage: install.R <yml_file>
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)

yml_file=args[1]

conda_path=conda_binary() #sometimes r environment doesn't recognize the conda path

# Specify the Conda environment you want to check
env_name <- "GenomeProt_env"


# Run the Conda command to list environments and capture the output
env_list <- system(paste0(conda_path," env list"),intern = TRUE)

#library path




# Check if the environment is in the list
if (any(grepl(env_name, env_list))) {
  message("The environment '", env_name, "' is present.")
  
} else {
  message("The environment '", env_name, "' is not found.")
  
   create_env_command <- paste("conda env create -f ",yml_file)
   print(create_env_command)
   system(create_env_command)
   
   #library path
   env_list_new <- system(paste0(conda_path," env list"),intern = TRUE)
   env_name_path<-env_list_new[grep("GenomeProt_env", env_list_new)]
   env_name_path_lib=paste0(strsplit(env_name_path,"\\s+")[[1]][2],"/lib/R/library/")
   
  
   # Install the required R packages into the specific conda environment (version 1.7.8.1)
   message("Installing xgboost from cloud.r-project.org")
   command_xgboost_install <- paste0(conda_path, " run -n ", env_name, " Rscript -e ", shQuote(paste0("install.packages('xgboost', repos='https://cloud.r-project.org', lib='", env_name_path_lib, "')")))
   
   system(command_xgboost_install)
   
   message("Installing ggtranscript from GitHub")
   
   command_ggtranscript_install <- paste(conda_path, "run -n", env_name, "Rscript -e",shQuote(paste0("withr::with_libpaths('", env_name_path_lib,"', devtools::install_github('dzhang32/ggtranscript', force=TRUE))")))
  
   system(command_ggtranscript_install)
}

