# usage: Rscript reformat_peptide_results.R -d <peptide results directory> -s <proteomics search tool> -o <path to write the reformatted peptide results file>

reformat_spectronaut_data <- function(peptide_data, dataset_id) {
  peptide_data_df <- peptide_data %>%
                       dplyr::select(dplyr::contains("PG."), dplyr::contains("PEP.AllOccurringProteinAccessions"), dplyr::contains("EG.PrecursorId"))

  if (length(peptide_data_df) == 0) {
    return(NULL)
  }

  metadata <- peptide_data_df %>%
                dplyr::mutate(Peptide = sapply(stringr::str_split(stringr::str_replace_all(EG.PrecursorId, "\\[.*?\\]", ""), "_"), `[`, 2)) %>%
                dplyr::select(Peptide, PEP.AllOccurringProteinAccessions) %>%
                base::unique() %>%
                dplyr::mutate(Dataset_id = dataset_id)
  return(metadata)
}

reformat_fragpipe_data <- function(peptide_data, dataset_id) {
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

reformat_fragpipe_quant_data <- function(peptide_data, dataset_id) {
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

reformat_dia_nn_data <- function(peptide_data, dataset_id) {
  peptide_data_colnames <- colnames(peptide_data)
  if (!("Stripped.Sequence" %in% peptide_data_colnames) | !("Protein.Ids" %in% peptide_data_colnames)) {
    return(NULL)
  }

  peptide_data_flt <- peptide_data %>%
                        dplyr::select(Stripped.Sequence, Protein.Ids) %>%
                        dplyr::mutate(Peptide = Stripped.Sequence) %>%
                        dplyr::select(Peptide, Protein.Ids) %>%
                        base::unique() %>%
                        dplyr::mutate(Dataset_id = dataset_id)

  return(peptide_data_flt)
}

reformat_maxquant_data <- function(peptide_data, dataset_id) {
  peptide_data_colnames <- colnames(peptide_data)
  if (!("Sequence" %in% peptide_data_colnames) | !("Leading razor protein" %in% peptide_data_colnames) | !("Proteins" %in% peptide_data_colnames)) {
    return(NULL)
  }

  peptide_data_flt <- peptide_data %>%
                        dplyr::select(Sequence, `Leading razor protein`, Proteins) %>%
                        dplyr::rename(Peptide = Sequence) %>%
                        dplyr::rename(Protein = `Leading razor protein`) %>%
                        dplyr::rename(Mapped_proteins = Proteins) %>%
                        dplyr::mutate(
                          Mapped_proteins = stringr::str_split(Mapped_proteins, ";"),
                          Mapped_proteins = lapply(Mapped_proteins, function (x) {
                              if (is.null(x) || length(x) == 0 || all(is.na(x))) return("")
                              x
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
              help = "peptide results directory; must contain .txt, .csv or .tsv files", metavar = "character"),
  make_option(c("-s", "--search_tool"), type = "character", default = NULL,
              help = "search tool used to generate the proteomic data (case-insensitive): Spectronaut / FragPipe (peptide.tsv) / FragPipe_quant (report.pr_matrix.tsv from FragPipe) / DIA-NN (report.pr_matrix.tsv from DIA-NN) / MaxQuant", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "path to write the reformatted peptide results file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# define inputs
peptide_results_dir <- opt$data_directory
search_tool <- opt$search_tool
output_path <- opt$output

# check if any of the arguments are missing

error_msg <- ""

if (is.null(peptide_results_dir)) {
  error_msg <- "Please provide a peptide results directory."
} else if (is.null(search_tool)) {
  error_msg <- "Please provide the search tool used to generate the proteomic data."
} else if (is.null(output_path)) {
  error_msg <- "Please provide a path to write the reformatted peptide results file to."
}

if (nchar(error_msg) > 0) {
  stop(error_msg)
}

# check the input arguments

search_tool <- tolower(trimws(as.character(search_tool)))
output_dir <- dirname(output_path)

if (!file.exists(peptide_results_dir)) {
  error_msg <- paste0("The peptide results directory '", peptide_results_dir, "' does not exist.")
} else if (!dir.exists(peptide_results_dir)) {
  error_msg <- paste0("The peptide results directory '", peptide_results_dir, "' exists but is not a directory.")
} else if (!(search_tool %in% c("spectronaut", "fragpipe", "fragpipe_quant", "dia-nn", "maxquant"))) {
  error_msg <- paste0("The search tool '", opt$search_tool, "' is not supported. Please specify one of the following (case-insensitive): 'Spectronaut', 'FragPipe' (peptide.tsv), 'FragPipe_quant' (report.pr_matrix.tsv from FragPipe), 'DIA-NN' (report.pr_matrix.tsv from DIA-NN), 'MaxQuant'")
} else if (file.exists(output_dir) && !dir.exists(output_dir)) {
  error_msg <- paste0("The directory of the output path '", output_dir, "' exists but is not a directory.")
} else {
  peptide_files_list <- list.files(path = peptide_results_dir, pattern = "\\.(txt|csv|tsv)$", all.files = TRUE, full.names = TRUE)
  peptide_files <- c()
  for (fn in peptide_files_list) {
    if (file.exists(fn) & (!dir.exists(fn)) & (file.size(fn) > 0)) {
      peptide_files <- c(peptide_files, fn)
    }
  }

  if (length(peptide_files) == 0) {
    error_msg <- paste0("There are no non-empty .txt, .csv or .tsv data files in the peptide results directory '", peptide_results_dir, "'.")
  }
}

if (nchar(error_msg) > 0) {
  stop(error_msg)
}

dataset_ids <- stringr::str_replace(basename(peptide_files), "\\.(txt|csv|tsv)$", "")
if (length(unique(dataset_ids)) != length(dataset_ids)) {
  stop(paste0("There are duplicate dataset IDs found in the peptide results directory '", peptide_results_dir, "'. Please ensure all dataset IDs are unique."))
} else if (any(nchar(dataset_ids) == 0)) {
  stop(paste0("Empty dataset IDs found in the peptide results directory '", peptide_results_dir, "'. Please ensure all dataset IDs consist of at least one character."))
}

# create the output directory if it does not already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
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
  dataset_id <- stringr::str_replace(basename(fn), "\\.(txt|csv|tsv)$", "")

  # reformat peptide results according to the proteomics search tool specified
  peptide_data <- readr::read_delim(fn)
  if (search_tool == "spectronaut") {
    results <- reformat_spectronaut_data(peptide_data, dataset_id)
  } else if (search_tool == "fragpipe") {
    results <- reformat_fragpipe_data(peptide_data, dataset_id)
  } else if (search_tool == "fragpipe_quant") {
    results <- reformat_fragpipe_quant_data(peptide_data, dataset_id)
  } else if (search_tool == "dia-nn") {
    results <- reformat_dia_nn_data(peptide_data, dataset_id)
  } else {
    results <- reformat_maxquant_data(peptide_data, dataset_id)
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
if (search_tool %in% c("spectronaut", "dia-nn")) {
  # pick whichever accession column is present
  acc_col <- intersect(c("PEP.AllOccurringProteinAccessions", "Protein.Ids"), colnames(combined_peptide_results))[1]
  if (is.na(acc_col)) {
    stop("Neither the 'PEP.AllOccurringProteinAccessions' column nor the 'Protein.Ids' column are found in the data.")
  }

  combined_peptide_results_uniq <- combined_peptide_results %>%
                                     dplyr::rename(Accessions = dplyr::all_of(acc_col)) %>%
                                     tidyr::separate_rows(Accessions, sep = ";") %>%
                                     dplyr::mutate(Accessions = trimws(Accessions)) %>%
                                     dplyr::filter(!is.na(Accessions), Accessions != "") %>%
                                     dplyr::group_by(Peptide) %>%
                                     dplyr::summarise(
                                       Accessions = paste(base::unique(trimws(Accessions)), collapse = ","),
                                       Dataset_id = paste(base::unique(trimws(Dataset_id)), collapse = ","),
                                       .groups = "drop"
                                     )

  uniq_peptides <- combined_peptide_results_uniq %>%
                     dplyr::mutate(
                       Protein = stringr::word(Accessions, 1, sep = stringr::fixed(",")),
                       Mapped_proteins = stringr::str_remove(Accessions, "^[^,]*,?")
                     ) %>%
                     dplyr::select(Peptide, Protein, Mapped_proteins, Dataset_id) %>%
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

readr::write_tsv(uniq_peptides, output_path)
message(paste0("Peptide results reformatting complete. Reformatted data written to '", output_path, "'."))