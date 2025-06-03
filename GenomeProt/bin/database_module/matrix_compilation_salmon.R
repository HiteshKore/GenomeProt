library(tximport)
library(optparse)
library(rtracklayer)
library(dplyr)
library(readr)



option_list = list(
  make_option(c("-q", "--salmon_outdir"), type="character", default=NULL, 
              help="Salmon output directory", metavar="character"),
  make_option(c("-s", "--sample_ids"), type="character", default=NULL, 
              help=" vector with sample ids", metavar="character"),
  make_option(c("-g", "--reference_gtf"), type="character", default=NULL, 
              help=" vector with sample ids", metavar="character")
); 

 opt_parser = OptionParser(option_list=option_list);
 opt = parse_args(opt_parser);

#print(opt$salmon_outdir)
samples=unlist(strsplit(opt$sample_ids, " "))

samples <- samples[samples != ""]

# # set path to salmon quant files
files <- file.path(opt$salmon_outdir, samples, "quant.sf")
print(files)
names(files) <- samples

# check if the files exist
print(all(file.exists(files)))

# use tximport to import salmon quantification files

txi <- tximport(files, type = "salmon", txOut = TRUE)


# import gtf to add gene information
gtf_data <- rtracklayer::import(opt$reference_gtf, format = "gtf")

# convert GTF data to a data frame
gtf_df <- as.data.frame(gtf_data)

# filter for relevant columns
transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
colnames(transcript_gene_info) <- c("TXNAME","GENEID")

# convert counts object to df
count_df <- as.data.frame(txi$counts)

# set col names
count_df <-  count_df %>% mutate(TXNAME = rownames(count_df)) %>% dplyr::select(TXNAME, everything())
rownames(count_df) <- NULL
# merge tx counts and gene info
count_df_merged <- left_join(count_df, transcript_gene_info, by="TXNAME")
count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, everything())

# export short-read count data
write_tsv(count_df_merged, file = paste0(opt$salmon_outdir, "counts_transcript.txt"), escape = "none", col_names = TRUE)