# usage: Rscript run_bambu.R -b <BAM file directory> -g <transcript annotation file (GTF/GFF) -s <organism> -o <output directory> -t <num_cpu_threads>

run_bambu_function <- function(bam_file_list, gtf, genomedb, output_directory, threads) {
  bambuAnnotations <- bambu::prepareAnnotations(gtf)
  se <- bambu::bambu(reads = bam_file_list, annotations = bambuAnnotations, genome = genomedb, verbose = TRUE, ncore = threads)
  bambu::writeBambuOutput(se, path = output_directory)

  tx_data <- as.data.frame(mcols(se))
  tx_data <- as.data.frame(apply(tx_data, 2, as.character))
  tx_data <- tx_data %>% dplyr::filter(novelTranscript == "TRUE")

  write.csv(tx_data, file.path(output_directory, "novel_transcript_classes.csv"), row.names = F, quote = F)
}

suppressPackageStartupMessages({
  library(optparse)
})

# define options
option_list = list(
  make_option(c("-b", "--bam"), type = "character", default = NULL,
              help = "BAM file directory", metavar = "character"),
  make_option(c("-g", "--gtf"), type = "character", default = NULL,
              help = "Transcript annotation file (GTF/GFF)", metavar = "character"),
  make_option(c("-s", "--species"), type = "character", default = NULL,
              help = "Organism: HUMAN / CAEEL / DROME / MOUSE / RAT / DANRE", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory", metavar = "character"),
  make_option(c("-t", "--thread"), type = "character", default = NULL,
              help = "Number of CPU threads to run bambu with", metavar = "integer")
);

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# store the input arguments
bam_dir <- opt$bam
gtf_file <- opt$gtf
organism <- opt$species
outdir <- opt$outdir
num_threads <- opt$thread

# check if any input arguments are missing

error_message <- ""

if (is.null(bam_dir)) {
  error_message <- "Please specify a directory containing BAM files."
} else if (is.null(gtf_file)) {
  error_message <- "Please provide a transcript annotation file (GTF/GFF)."
} else if (is.null(organism)) {
  error_message <- "Please specify one of the following organisms (case-insensitive): 'HUMAN' (H. sapiens), 'CAEEL' (C. elegans), 'DROME' (D. melanogaster), 'MOUSE' (M. musculus), 'RAT' (R. norvegicus) or 'DANRE' (D. rerio)."
} else if (is.null(outdir)) {
  error_message <- "Please specify an output directory."
} else if (is.null(num_threads)) {
  error_message <- "Please specify the number of CPU threads to run bambu with."
}

if (nchar(error_message) > 0) {
  stop(error_message)
}

# check the input arguments

organism <- toupper(trimws(as.character(organism)))
num_threads <- floor(as.numeric(num_threads))
bam_files <- NULL

if (!file.exists(bam_dir)) {
  error_message <- paste0("The BAM file directory '", bam_dir, "' does not exist.")
} else if (file.exists(bam_dir) && !dir.exists(bam_dir)) {
  error_message <- paste0("'", bam_dir, "' exists but is not a directory.")
} else if (!file.exists(gtf_file)) {
  error_message <- paste0("The transcript annotation file '", gtf_file, "' does not exist.")
} else if (!(organism %in% c("HUMAN", "CAEEL", "DROME", "MOUSE", "RAT", "DANRE"))) {
  error_message <- paste0("Organism '", as.character(opt$species), "' is not supported. Please specify one of the following (case-insensitive): 'HUMAN' (H. sapiens), 'CAEEL' (C. elegans), 'DROME' (D. melanogaster), 'MOUSE' (M. musculus), 'RAT' (R. norvegicus) or 'DANRE' (D. rerio).")
} else if (file.exists(outdir) && !dir.exists(outdir)) {
  error_message <- paste0("'", outdir, "' exists but is not a directory.")
} else if (!(is.finite(num_threads) && num_threads > 0)) {
  error_message <- paste0("The number of CPU threads '", opt$thread, "' cannot be parsed as a positive integer.")
} else {
  bam_files <- list.files(path = bam_dir, "\\.bam$", full.names = TRUE)
  if (length(bam_files) == 0) {
    error_message <- paste0("The BAM file directory '", bam_dir, "' does not contain any BAM files (i.e. files ending with '.bam').")
  }
}

if (nchar(error_message) > 0) {
  stop(error_message)
}

# create the output directory if it does not already exist
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# load the rest of the libraries needed for the remainder of the script to run
suppressPackageStartupMessages({
  library(bambu)
  library(dplyr)
  library(Rsamtools)
})

# create list of BAMs
bam_file_list <- Rsamtools::BamFileList(bam_files)

# remove the .bam file extension
bam_file_names <- list.files(path = bam_dir, "\\.bam$", full.names = FALSE)

# rename the list of BAMs to their original names
names(bam_file_list) <- sub("\\.bam$", "", bam_file_names)
print(bam_file_list)

# select the reference genome corresponding to the specified organism
if (organism == "HUMAN") {          # H. sapiens
  suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
  })
  genomedb <- BSgenome.Hsapiens.UCSC.hg38
} else if (organism == "CAEEL") {   # C. elegans
  suppressPackageStartupMessages({
    library(BSgenome.Celegans.UCSC.ce11)
  })
  genomedb <- BSgenome.Celegans.UCSC.ce11
} else if (organism == "DROME") {   # D. melanogaster
  suppressPackageStartupMessages({
    library(BSgenome.Dmelanogaster.UCSC.dm6)
  })
  genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
} else if (organism == "MOUSE") {   # M. musculus
  suppressPackageStartupMessages({
    library(BSgenome.Mmusculus.UCSC.mm39)
  })
  genomedb <- BSgenome.Mmusculus.UCSC.mm39
} else if (organism == "RAT") {     # R. norvegicus
  suppressPackageStartupMessages({
    library(BSgenome.Rnorvegicus.UCSC.rn7)
  })
  genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
} else {                            # D. rerio
  suppressPackageStartupMessages({
    library(BSgenome.Drerio.UCSC.danRer11)
  })
  genomedb <- BSgenome.Drerio.UCSC.danRer11
}
#} else if (organism == "DANRE") {   # D. rerio
#  library(BSgenome.Drerio.UCSC.danRer11)
#  genomedb <- BSgenome.Drerio.UCSC.danRer11
#} else if (organism == "PANTR") {   # P. troglodytes
#  library(BSgenome.Ptroglodytes.UCSC.panTro6)
#  genomedb <- BSgenome.Ptroglodytes.UCSC.panTro6
#} else if (organism == "BOVIN") {   # B. taurus
#  library(BSgenome.Btaurus.UCSC.bosTau9)
#  genomedb <- BSgenome.Btaurus.UCSC.bosTau9
#} else if (organism == "XENTR") {   # X. tropicalis
#  library(...)
#  genomedb <- ...
#} else {                            # S. cerevisiae
#  library(BSgenome.Scerevisiae.UCSC.sacCer3)
#  genomedb <- BSgenome.Scerevisiae.UCSC.sacCer3
#}

# run bambu
run_bambu_function(bam_file_list, gtf_file, genomedb, outdir, num_threads)