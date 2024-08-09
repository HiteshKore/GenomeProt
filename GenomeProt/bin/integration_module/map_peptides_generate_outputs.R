
# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

suppressMessages({
  library(data.table)
  library(ORFik)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringdist)
  library(Repitools)
  library(optparse)
  library(mygene)
})

source("global.R")
source("R/integration_functions.R")

option_list = list(
  make_option(c("-p", "--proteomics"), type="character", default=NULL,
              help="Proteomics data file", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Custom FASTA used for proteomics", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF used to generate custom FASTA", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

system("mkdir integ_output")

proteomics_import_file <- opt$proteomics
fasta_import_file <- opt$fasta
gtf_import_file <- opt$gtf

# source("~/Documents/GenomeProt_tmp/GenomeProt/GenomeProt/R/integration_functions.R")
# proteomics_import_file <- "~/Documents/linda_data/genomeprot_peptides.txt"
# fasta_import_file <- "~/Documents/linda_data/proteome_database.fasta"
# gtf_import_file <- "~/Documents/linda_data/proteome_database_transcripts.gtf"

# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))

gtf <- makeTxDbFromGFF(gtf_import_file) # make txdb of gtf

md <- import_fasta(fasta_import_file, pd, gtf_import_file)

# ---------------------------------------- #


# ------------- apply functions ------------- #

# extract ORF transcript coordinates to df
md$orf_tx_id <- paste0(md$protein_name, "_", md$transcript_id)

# get just unique orf and transcript for mapping
orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id, strand)
orf_transcript_coords_df <- orf_transcript_coords_df[!(base::duplicated(orf_transcript_coords_df)),]

# make GRanges from df of ORF transcript coordinates
orf_transcript_coords <- makeGRangesFromDataFrame(orf_transcript_coords_df,
  keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
  seqnames.field="transcript_id", start.field="txstart", end.field="txend", strand.field="strand",
  starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(orf_transcript_coords) <- c(orf_transcript_coords$orf_tx_id) # set names
mcols(orf_transcript_coords)$gene_id <- c(orf_transcript_coords_df$gene_id)

# get exons for mapping
exons <- exonsBy(gtf, "tx", use.names=TRUE) # get exon data per transcript
exons_filt <- exons[names(exons) %in% orf_transcript_coords_df$transcript_id] # filter for only transcripts with peptides

orf_tx_names <- as.character(seqnames(orf_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(orf_transcript_coords) <- match(orf_tx_names, names(exons_filt)) 

# original peptides in the GRanges
orf_ids <- orf_transcript_coords$orf_tx_id

# ORFik map to genome coordinates
orf_in_genomic <- ORFik::pmapFromTranscriptF(orf_transcript_coords, exons_filt, removeEmpty = T)

# map back to GRangesList, with group information
orf_in_genomic@unlistData$PID <- orf_ids[groupings(orf_in_genomic)]

# add exon_number for GTF export
orf_in_genomic_gr <- unlist(orf_in_genomic, use.names=F) # convert to GRanges

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(orf_in_genomic_gr), mcols(orf_in_genomic_gr)$PID, FUN = seq_along)
# add to GRanges
mcols(orf_in_genomic_gr)$exon_number <- exon_number_vec
# re-list
orf_in_genomic <- split(orf_in_genomic_gr, ~ mcols(orf_in_genomic_gr)$PID)

# peptides
# use ORF transcript coords to determine peptide transcript coords
peptide_transcript_coords <- extract_peptide_coords(md)

# map to genomic coords
peptide_tx_names <- as.character(seqnames(peptide_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(peptide_transcript_coords) <- match(peptide_tx_names, names(exons_filt)) 

# original peptides in the GRanges
pep_ids <- peptide_transcript_coords$peptide
pep_PID_ids <- peptide_transcript_coords$PID
pep_gene_ids <- peptide_transcript_coords$gene_id

# ORFik map to genome coordinates

# causes script to fail
pep_in_genomic <- ORFik::pmapFromTranscriptF(peptide_transcript_coords, exons_filt, removeEmpty = F)

# map back to GRangesList, with group information
pep_in_genomic@unlistData$peptide <- pep_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$PID <- pep_PID_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$gene <- pep_gene_ids[groupings(pep_in_genomic)]

# remove 0 ranges
pep_in_genomic_gr <- unlist(pep_in_genomic, use.names=F) # convert to GRanges
# rename with transcript and peptide
tx_pep_names <- c(paste0(names(pep_in_genomic_gr), "_", pep_in_genomic_gr$peptide))
names(pep_in_genomic_gr) <- tx_pep_names # set names
pep_in_genomic_gr$tx_pid_grouping <- paste0(pep_in_genomic_gr$PID, "_", names(pep_in_genomic_gr))
# remove 0 ranges
pep_in_genomic_gr <- subset(pep_in_genomic_gr, (start(pep_in_genomic_gr) != 0 & end(pep_in_genomic_gr) != 0)  )

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(pep_in_genomic_gr), pep_in_genomic_gr$tx_pid_grouping, FUN = seq_along)
# add to GRanges
mcols(pep_in_genomic_gr)$exon_number <- exon_number_vec
mcols(pep_in_genomic_gr)$tx_pid_grouping <- NULL

# re-list
pep_in_genomic <- split(pep_in_genomic_gr, ~ names(pep_in_genomic_gr))

# ---------------------------------------- #


# ------------- export ------------- #

# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, "integ_output/peptides.bed12", rgb = 0)

# export bed12 of ORFs
# should we change it so that only unique ORFs are exported, not every ORF per transcript?
ORFik::export.bed12(orf_in_genomic, "integ_output/ORFs.bed12", rgb = 0)

# format GTF of ORFs
orf_in_genomic_gr$source <- c("custom")
orf_in_genomic_gr$type <- c("CDS")
orf_in_genomic_gr$phase <- 0
orf_in_genomic_gr$ORF_id <- names(orf_in_genomic_gr)
orf_in_genomic_gr$transcript_id <- names(orf_in_genomic_gr)
names(orf_in_genomic_gr) <- NULL
orf_in_genomic_gr$group_id <- "ORFs"

# format GTF of all transcripts that had mapped peptides
gtf_for_exporting <- import(gtf_import_file, format="gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]
gtf_filtered$group_id <- "transcripts"

# reformat exons for bed12
gtf_as_bed12 <- gtf_filtered[mcols(gtf_filtered)$type == "exon"]

names(gtf_as_bed12) <- paste0(gtf_as_bed12$transcript_id, "_", gtf_as_bed12$gene_id)

# convert to grl
tx_in_genomic <- split(gtf_as_bed12, ~ names(gtf_as_bed12))

# export bed12 of transcripts
ORFik::export.bed12(tx_in_genomic, "integ_output/transcripts.bed12", rgb = 0)

# ---------------------------------------- #


# ------- summary file -------- #

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
results_pept_df <- pep_in_genomic_gr %>% as_tibble()
results_pept_df <- separate(results_pept_df, txname, into = c("transcript_id", "peptide"), sep = "_", remove = TRUE)
results_pept_df$gene_id <- results_pept_df$gene
results_pept_df$gene <- NULL

# group by peptide and transcript to summarise based on how many exons peptide spans
results_pept_df_unique <- results_pept_df %>% 
  dplyr::group_by(peptide, transcript_id) %>% 
  dplyr::slice_max(exon_number) %>% 
  dplyr::mutate(number_exons = exon_number) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end, -width, -exon_number)

# merge results with metadata
metadata_to_merge <- md %>% 
  dplyr::select(PID, peptide, transcript_id, gene_id, gene_name)

peptide_result <- merge(results_pept_df_unique, metadata_to_merge, by=c("PID", "peptide", "gene_id", "transcript_id"), all.x=T, all.y=F)

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

# define transcripts as known or novel
# currently, python script requires novel txs to have prefix 'Bambu'
peptide_result <- peptide_result %>% 
  mutate(isoform_status = case_when(
    #!startsWith(transcript_id, "EN") ~ "novel",
    startsWith(transcript_id, "Bambu") ~ "novel",
    TRUE ~ "known"
  ))

# determine ORF and peptide mapping status

# peptides mapped to a unique ORF (and gene) = high conf
# peptides mapped to multi ORFs from one gene = high conf (unique genomic location is identified)
# peptides mapped to multi genes locations = low conf

peptide_result <- peptide_result %>% 
  dplyr::group_by(peptide) %>% 
  dplyr::mutate(peptide_status = case_when(
    length(unique(gene_id))==1 ~ "high",
    length(unique(PID))>1 & length(unique(gene_id))>1 ~ "low",
    TRUE ~ NA),
    orf_status_1 = case_when(
      length(unique(gene_id))==1 & length(unique(PID))==1 ~ "high",
      TRUE ~ "low")) %>% 
  ungroup()

# if an ORF is identified by at least one high confidence peptide, then the ORF is high confidence
peptide_result <- peptide_result %>% 
  dplyr::group_by(PID) %>% 
  dplyr::mutate(orf_status = case_when(
    orf_status_1 %in% c("high") ~ "high",
    TRUE ~ "low")) %>% 
  dplyr::ungroup() %>% dplyr::select(-orf_status_1)

peptide_result <- peptide_result %>% 
  dplyr::group_by(peptide) %>% 
  dplyr::mutate(iso_map_status_1 = case_when(
    peptide_status %in% c("high") & length(unique(transcript_id))==1 & startsWith(transcript_id, "Bambu") ~ "novel_iso_uniq_identified",
    peptide_status %in% c("high") & length(unique(transcript_id))==1 & !startsWith(transcript_id, "Bambu") ~ "known_iso_uniq_identified",
    TRUE ~ NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(transcript_id) %>% 
  dplyr::mutate(iso_map_status = case_when(
    iso_map_status_1 %in% c("novel_iso_uniq_identified") ~ "novel_iso_uniq_identified",
    iso_map_status_1 %in% c("known_iso_uniq_identified") ~ "known_iso_uniq_identified",
    TRUE ~ "other")) %>% 
  dplyr::ungroup() %>% dplyr::select(-iso_map_status_1)

# include orf_status and peptide_status in GTF mcols
results_pept_df$peptide_status <- NULL
results_to_merge_with_granges <- merge(results_pept_df, peptide_result, by=c("transcript_id", "peptide", "strand", "PID", "gene_id", "seqnames"), all.x=T, all.y=F)
results_to_merge_with_granges <- results_to_merge_with_granges[!(duplicated(results_to_merge_with_granges)),]
results_to_merge_with_granges$naming <- paste0(results_to_merge_with_granges$transcript_id, "_", results_to_merge_with_granges$peptide)

# make GRanges from df of ORF transcript coordinates
pep_in_genomic_gr_export <- makeGRangesFromDataFrame(results_to_merge_with_granges,
                                                  keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                  seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand",
                                                  starts.in.df.are.0based=FALSE, na.rm=TRUE)
names(pep_in_genomic_gr_export) <- c(pep_in_genomic_gr_export$naming) # set names

# add mcols
pep_in_genomic_gr_export$source <- "custom"
pep_in_genomic_gr_export$type <- "exon"
pep_in_genomic_gr_export$group_id <- "peptides"

# export summary data
write.csv(peptide_result, "integ_output/peptide_info.csv", row.names=F, quote=F)

# export annotations for vis
combined <- c(pep_in_genomic_gr_export, orf_in_genomic_gr, gtf_filtered)
export(combined, "integ_output/combined_annotations.gtf", format="gtf")

# ---------------------------------------- #



