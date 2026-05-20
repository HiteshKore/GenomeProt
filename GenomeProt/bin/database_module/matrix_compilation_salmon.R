# usage: Rscript matrix_compilation_salmon.R -q <Salmon output directory> -s <list of space-separated sample IDs> -g <reference GTF file>

# define options
library(optparse)
option_list = list(
  make_option(c("-q", "--salmon_outdir"), type = "character", default = NULL,
              help = "Salmon output directory", metavar = "character"),
  make_option(c("-s", "--sample_ids"), type = "character", default = NULL,
              help = "List of space-separated sample IDs", metavar = "character"),
  make_option(c("-g", "--reference_gtf"), type = "character", default = NULL,
              help = "Reference GTF file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# store the input arguments
salmon_outdir <- opt$salmon_outdir
sample_ids <- opt$sample_ids
reference_gtf <- opt$reference_gtf

# check if any input arguments are missing
if (is.null(salmon_outdir)) {
  stop("Please specify a Salmon output directory.")
} else if (is.null(sample_ids)) {
  stop("Please provide a list of space-separated sample IDs.")
} else if (is.null(reference_gtf)) {
  stop("Please provide a reference GTF file.")
}

# check the input arguments

samples <- unlist(strsplit(sample_ids, " "))
samples <- samples[samples != ""]

if (length(samples) == 0) {
  stop("The list of sample IDs provided is empty or consists entirely of whitespace.")
} else if (!file.exists(salmon_outdir)) {
  stop(paste0("The Salmon output directory '", salmon_outdir, "' does not exist."))
} else if (file.exists(salmon_outdir) && !dir.exists(salmon_outdir)) {
  stop(paste0("'", salmon_outdir, "' exists but is not a directory."))
} else if (!file.exists(reference_gtf)) {
  stop(paste0("The reference GTF file '", reference_gtf, "' does not exist."))
}

library(tximport)
library(rtracklayer)
library(dplyr)
library(readr)

# set path to salmon quant files
files <- file.path(salmon_outdir, samples, "quant.sf")
names(files) <- samples

# use tximport to import salmon quantification files
txi <- tximport(files, type = "salmon", txOut = TRUE)

# import the reference GTF to add gene information
gtf_df <- as.data.frame(rtracklayer::import(reference_gtf, format = "gtf"))

# filter for relevant columns
transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
colnames(transcript_gene_info) <- c("TXNAME", "GENEID")

# convert counts object to df
count_df <- as.data.frame(txi$counts)

# set col names
count_df <- count_df %>% mutate(TXNAME = rownames(count_df)) %>% dplyr::select(TXNAME, everything())
rownames(count_df) <- NULL

# merge tx counts and gene info
count_df_merged <- left_join(count_df, transcript_gene_info, by = "TXNAME")
count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, everything())

# export short-read count data
write_tsv(count_df_merged, file = file.path(salmon_outdir, "counts_transcript.txt"), escape = "none", col_names = TRUE)