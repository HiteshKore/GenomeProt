# usage: Rscript generate_integration_summary_report.R -i <integration summary report Rmd file> -o <output directory>

# define options
library(optparse)
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "integration summary report Rmd file", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# store the input arguments
report_rmd_file <- opt$input
report_outdir <- opt$outdir

# check the input arguments
if (!file.exists(report_rmd_file)) {
  stop(paste0("The integration summary report Rmd file '", report_rmd_file, "' does not exist."))
} else if (!file.exists(report_outdir)) {
  stop(paste0("The output directory '", report_outdir, "' does not exist."))
} else if (file.exists(report_outdir) && !dir.exists(report_outdir)) {
  stop(paste0("'", report_outdir, "' exists but is not a directory."))
} else if (!file.exists(file.path(report_outdir, "peptide_info.tsv"))) {
  stop(paste0("The file 'peptide_info.tsv' cannot be found in the output directory '", report_outdir, "'."))
}

# create the directory for storing the report images if it does not exist
report_outdir_image <- file.path(report_outdir, "report_images")
if (!dir.exists(report_outdir_image)) {
  dir.create(report_outdir_image)
}

# generate the report
library(rmarkdown)
rmarkdown::render(input = report_rmd_file,
                  output_file = file.path(report_outdir, "summary_report.html"),
                  output_format = "html_document",
                  params = list(directory = report_outdir, file = "peptide_info.tsv"),
                  intermediates_dir = report_outdir)