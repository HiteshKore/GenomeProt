#!/bin/bash
# usage: install.R <yml_file>
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)

yml_file = args[1]

conda_path = conda_binary() # get the path to the conda binary

# Specify the Conda environment you want to check
env_name <- "GenomeProt_env"

# Run the Conda command to list environments and capture the output
env_list <- system(paste0(conda_path," env list"), intern = TRUE)

# Check if the environment is in the list
if (any(grepl(env_name, env_list))) {
  message("The environment '", env_name, "' is present.")
} else {
  message("The environment '", env_name, "' is not found.")

  create_env_command <- paste("conda env create -f ",yml_file)
  print(create_env_command)
  system(create_env_command)
}