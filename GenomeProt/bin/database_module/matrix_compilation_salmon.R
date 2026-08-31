# usage: Rscript matrix_compilation_salmon.R -s <Salmon output directory> -g <reference GTF file> -o <path to write the counts file>

# define options
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-s", "--salmon_outdir"), type = "character", default = NULL,
              help = "Salmon output directory", metavar = "character"),
  make_option(c("-g", "--reference_gtf"), type = "character", default = NULL,
              help = "Reference GTF file", metavar = "character"),
  make_option(c("-o", "--counts_output_file"), type = "character", default = NULL,
              help = "Path to write the counts file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# store the input arguments
salmon_outdir <- opt$salmon_outdir
reference_gtf <- opt$reference_gtf
counts_output_file <- opt$counts_output_file

# check if any input arguments are missing
if (is.null(salmon_outdir)) {
  stop("Please specify a Salmon output directory.")
} else if (is.null(reference_gtf)) {
  stop("Please provide a reference GTF file.")
} else if (is.null(counts_output_file)) {
  stop("Please provide a path to write the counts file.")
}

salmon_outdir=list.dirs(salmon_outdir, full.names = TRUE, recursive = FALSE)

salmon_outdir <- salmon_outdir[basename(salmon_outdir) != "salmon_index"]

# check the input arguments
files <- Sys.glob(file.path(salmon_outdir,"/quant.sf"))

if (length(files) == 0) {
    stop(paste0("The Salmon output directory '", salmon_outdir, "' does not contain any directories with the file 'quant.sf'."))
}


if (all(file.exists(files))) {
  # load the rest of the libraries needed for the remainder of the script to run
  suppressPackageStartupMessages({
    library(tximport)
    library(rtracklayer)
    library(dplyr)
    library(readr)
  })
  
  # import and filter the GTF for relevant gene info
  gtf_df <- as.data.frame(rtracklayer::import(reference_gtf, format = "gtf"))
  transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
  colnames(transcript_gene_info) <- c("TXNAME", "GENEID")
  
  # import salmon quantification files
  
  sample <- basename(dirname(files))
  names(files) <- sample
  txi <- tximport::tximport(files, type = "salmon", txOut = TRUE)
  
  # prepare the transcript counts data
  count_df <- as.data.frame(txi$counts) %>% dplyr::mutate(TXNAME = rownames(txi$counts)) %>% dplyr::select(TXNAME, dplyr::everything())
  
  rownames(count_df) <- NULL
  
  # merge the transcript counts and gene info
  count_df_merged <- dplyr::left_join(count_df, transcript_gene_info, by = "TXNAME")
  count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, dplyr::everything())
  
  # export short-read count data
  readr::write_tsv(count_df_merged, file = counts_output_file, escape = "none", col_names = TRUE)

}else{
  
  stop(paste0("'quant.sf' is not present for some samples."))
}

