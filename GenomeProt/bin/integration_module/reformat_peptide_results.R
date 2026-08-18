# usage: Rscript reformat_peptide_results.R -d <peptide results directory> -s <proteomics search tool: one of 'Spectronaut', 'FragPipe' (peptide.tsv) or 'FragPipe_quant' (report.pr_matrix.tsv)>

reformat_spectronaut_data <- function(peptide_file, dataset_id) {
  peptide_data <- readr::read_tsv(peptide_file)
  peptide_data_df <- peptide_data %>%
                       dplyr::select(dplyr::contains("PG."), dplyr::contains("PEP.AllOccurringProteinAccessions"), dplyr::contains("EG.PrecursorId"))

  if (length(peptide_data_df) == 0) {
    return(NULL)
  }

  metadata <- peptide_data_df %>%
                dplyr::mutate(Peptide = sapply(stringr::str_split(stringr::str_replace_all(EG.PrecursorId, "\\[.*?\\]", ""), "_"), `[`, 2))
                dplyr::select(Peptide, PEP.AllOccurringProteinAccessions) %>%
                base::unique() %>%
                dplyr::mutate(Dataset_id = dataset_id)
  return(metadata)
}

reformat_fragpipe_data <- function(peptide_file, dataset_id) {
  peptide_data <- readr::read_tsv(peptide_file)

  peptide_data_colnames <- colnames(peptide_data)
  if (!("Peptide" %in% peptide_data_colnames) | !("Protein" %in% peptide_data_colnames) | !("Mapped Proteins" %in% peptide_data_colnames)) {
    return(NULL)
  }

  peptide_data_flt <- peptide_data %>%
                        dplyr::select(Peptide, Protein, `Mapped Proteins`) %>%
                        dplyr::mutate(
                          # Extract accession from Protein (middle piece between | |)
                          Protein = stringr::str_extract(Protein, "(?<=\\|)[^|]+(?=\\|)"),
                          `Mapped Proteins` = stringr::str_split(`Mapped Proteins`, ",\\s*"),
                          `Mapped Proteins` = lapply(`Mapped Proteins`, function(x) {
                              if (all(is.na(x))) return("")
                              stringr::str_extract(x, "(?<=\\|)[^|]+(?=\\|)")
                          }),
                          # Collapse vector back into comma-separated string
                          `Mapped Proteins` = sapply(`Mapped Proteins`, function(x) {
                              if (all(is.na(x))) return("")
                              paste(x, collapse = ",")
                          })
                        ) %>%
                        dplyr::rename(Mapped_proteins = `Mapped Proteins`) %>%
                        dplyr::mutate(Dataset_id = dataset_id)
  return(peptide_data_flt)
}

reformat_fragpipe_quant_data <- function(peptide_file, dataset_id) {
  peptide_data <- readr::read_tsv(peptide_file)

  peptide_data_colnames <- colnames(peptide_data)
  if (!("Stripped.Sequence" %in% peptide_data_colnames) | !("Protein.Ids" %in% peptide_data_colnames) | !("All Mapped Proteins" %in% peptide_data_colnames)) { 
    return(NULL)       
  }

  peptide_data_flt <- peptide_data %>%
                        dplyr::select(Stripped.Sequence, Protein.Ids, `All Mapped Proteins`)
  colnames(peptide_data_flt) <- c("Peptide", "Protein", "Mapped_proteins")

  peptide_data_flt <- peptide_data_flt %>%
                        dplyr::mutate(
                          # Extract accession from Protein (middle piece between | |)
                          Protein = stringr::str_extract(Protein, "(?<=\\|)[^|]+(?=\\|)"),
                          Mapped_proteins = stringr::str_split(Mapped_proteins, ","),
                          Mapped_proteins = lapply(Mapped_proteins, function(x) {
                              if (all(is.na(x))) return("")
                              stringr::str_extract(x, "(?<=\\|)[^|]+(?=\\|)")
                          }),
                          Mapped_proteins = mapply(function(mapped, prot) {
                              if (is.null(mapped) || length(mapped) == 0) return("")
                              mapped <- mapped[!is.na(mapped)]  # remove NAs
                              mapped <- trimws(mapped)          # remove spaces
                              mapped <- mapped[mapped != prot]  # remove the protein accession
                              if (length(mapped) == 0) "" else paste(mapped, collapse = ",")
                          }, Mapped_proteins, Protein),
                          Mapped_proteins = as.character(Mapped_proteins),
                          Dataset_id = dataset_id
                        )
  return(peptide_data_flt)
}

suppressPackageStartupMessages({
  library(optparse)
})

# define options
option_list = list(
  make_option(c("-d", "--data_directory"), type = "character", default = NULL,
              help = "peptide results directory; must contain .tsv files (NOTE: peptide_data.tsv will be renamed to peptide_data_renamed.tsv as reformatted data is written to peptide_data.tsv)", metavar = "character"),
  make_option(c("-s", "--search_tool"), type = "character", default = NULL,
              help = "search tool used to generate the proteomic data: specify 'Spectronaut', 'FragPipe' (peptide.tsv), or 'FragPipe_quant' (report.pr_matrix.tsv)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# define inputs
peptide_results_dir <- opt$data_directory
search_tool <- opt$search_tool

# check if any of the arguments are missing

error_msg <- ""

if (is.null(peptide_results_dir)) {
  error_msg <- "Please provide a peptide results directory."
} else if (is.null(search_tool)) {
  error_msg <- "Please provide the search tool used to generate the proteomic data."
}

if (nchar(error_msg) > 0) {
  stop(error_msg)
}

# check the input arguments

search_tool <- tolower(trimws(as.character(search_tool)))
if (!file.exists(peptide_results_dir)) {
  error_msg <- paste0("The peptide results directory '", peptide_results_dir, "' does not exist.")
} else if (!dir.exists(peptide_results_dir)) {
  error_msg <- paste0("The peptide results directory '", peptide_results_dir, "' exists but is not a directory.")
} else if (!(search_tool %in% c("spectronaut", "fragpipe", "fragpipe_quant"))) {
  error_msg <- paste0("The search tool '", opt$search_tool, "' is not supported. Please specify one of the following (case-insensitive): 'Spectronaut', 'FragPipe' (peptide.tsv), 'FragPipe_quant' (report.pr_matrix.tsv)")
} else {
  peptide_files_list <- list.files(path = peptide_results_dir, pattern = "\\.tsv$", all.files = TRUE, full.names = TRUE)
  peptide_files <- c()
  for (fn in peptide_files_list) {
    if (file.exists(fn) & (!dir.exists(fn)) & (file.size(fn) > 0)) {
      peptide_files <- c(peptide_files, fn)
    }
  }

  if (length(peptide_files) == 0) {
    error_msg <- paste0("There are no non-empty .tsv data files in the peptide results directory '", peptide_results_dir, "'.")
  }
}

if (nchar(error_msg) > 0) {
  stop(error_msg)
}

# load the rest of the libraries needed for the remainder of the script to run
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(stats)
})

# process each peptide results file
peptide_results <- list()
for (fn in peptide_files) {
  dataset_id <- stringr::str_replace(basename(fn), "\\.tsv$", "")

  # rename peptide_data.tsv to peptide_data_renamed.tsv
  if (dataset_id == "peptide_data") {
    new_filename <- file.path(peptide_results_dir, "peptide_data_renamed.tsv")
    message(paste0("Renaming '", fn, "' to '", new_filename, "' before reformatting..."))
    file.rename(fn, new_filename)
    fn <- new_filename
    dataset_id <- "peptide_data_renamed"
  }

  # reformat peptide results according to the proteomics search tool specified
  if (search_tool == "spectronaut") {
    results <- reformat_spectronaut_data(fn, dataset_id)
  } else if (search_tool == "fragpipe") {
    results <- reformat_fragpipe_data(fn, dataset_id)
  } else {
    results <- reformat_fragpipe_quant_data(fn, dataset_id)
  }

  if (length(results) > 0) {
    message(paste0("Successfully parsed '", fn, "'."))
    peptide_results[[dataset_id]] <- results
  } else {
    message(paste0("Warning: Unable to parse '", fn, "' due to missing column names. Check if the file is formatted correctly and if the correct search tool is provided. Skipping..."))
  }
}

combined_peptide_results <- do.call(rbind, peptide_results)
if (length(combined_peptide_results) == 0) {
  error_msg <- "Unable to parse any file in the peptide results directory. Please ensure that the files are formatted correctly and the correct search tool is provided. Stopping..."
  stop(error_msg)
}

# combine peptide results
if (search_tool == "spectronaut") {
  combined_peptide_results_uniq <- combined_peptide_results %>%
                                     tidyr::separate_rows(PEP.AllOccurringProteinAccessions, sep = ";") %>%
                                     dplyr::group_by(Peptide) %>%
                                     dplyr::summarise(
                                       PEP.AllOccurringProteinAccessions = paste(base::unique(trimws(PEP.AllOccurringProteinAccessions)), collapse = ","),
                                       Dataset_id = paste(base::unique(trimws(Dataset_id)), collapse = ","),
                                       .groups = "drop"
                                     )

  uniq_peptides <- combined_peptide_results_uniq %>%
    dplyr::mutate(Protein = sapply(stringr::str_split(PEP.AllOccurringProteinAccessions, ","), `[`, 1)) %>%
    dplyr::rename(Mapped_proteins = PEP.AllOccurringProteinAccessions) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Mapped_proteins = paste(stringr::str_split(Mapped_proteins, ",")[[1]][-1], collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::select("Peptide", "Protein", "Mapped_proteins", "Dataset_id") %>%
    stats::na.omit()
} else {
  uniq_peptides <- combined_peptide_results %>%
    tidyr::separate_rows(Mapped_proteins, sep = ",") %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(
      Protein = paste(base::unique(trimws(Protein)), collapse = ","),
      Mapped_proteins = paste(base::unique(trimws(Mapped_proteins)), collapse = ","),
      Dataset_id = paste(base::unique(trimws(Dataset_id)), collapse = ","),
      .groups = "drop"
    )
}

peptide_data_filename <- file.path(peptide_results_dir, "peptide_data.tsv")
readr::write_tsv(uniq_peptides, peptide_data_filename)
message(paste0("Peptide results reformatting complete. Reformatted data written to '", peptide_data_filename, "'."))