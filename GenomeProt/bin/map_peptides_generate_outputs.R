

# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

# required columns in peptides.tsv data
# c(`Mapped Proteins`, Protein, Peptide)

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

#source("~/Documents/GenomeProt_tmp/GenomeProt/GenomeProt/R/functions.R")
source("R/functions.R")
source("global.R")

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

proteomics_import_file <- opt$proteomics
fasta_import_file <- opt$fasta
gtf_import_file <- opt$gtf

# proteomics_import_file <- "~/Documents/pg_server/integ_test_data/report.pr_matrix.tsv"
# fasta_import_file <- "~/Documents/pg_server/integ_test_data/ProteomeDb.fasta"
# gtf_import_file <- "~/Documents/pg_server/integ_test_data/ORFome_transcripts.gtf"

# ------------- args input ------------- #
# # replace with args[]
# #proteomics_import_file <- "~/Documents/proteogenomics/2024/miguel_tx_and_proteomics/Fragpipe_results/peptide.tsv"
# proteomics_import_file <- "~/Documents/proteogenomics/2024/ben_data/report.pr_matrix.tsv"
# 
# # change requirement to fasta
# metadata_import_file <- "~/Documents/proteogenomics/2024/ben_data/ProteomeDb_metadata.txt"
# fasta_import_file <- "~/Documents/proteogenomics/2024/ben_data/ProteomeDb.fasta"
# # still require GTF
# gtf_import_file <- "~/Documents/proteogenomics/2024/ben_data/ORFome_transcripts.gtf"
# 
# output_path <- "~/Documents/proteogenomics/2024/ben_data/peptide_results"

# -------------------------------------- #


# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))

gtf <- makeTxDbFromGFF(gtf_import_file) # make txdb of gtf

md <- import_fasta(fasta_import_file, pd, gtf)


# ---------------------------------------- #


# ------------- apply functions ------------- #

# extract ORF transcript coordinates to df
md$orf_tx_id <- paste0(md$protein_name, "_", md$transcript_id)

orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id)
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
#orf_gene_ids <- orf_transcript_coords$gene

# ORFik map to genome coordinates
orf_in_genomic <- ORFik::pmapFromTranscriptF(orf_transcript_coords, exons_filt, removeEmpty = T)

#orf_in_genomic <- mapFromTranscripts(orf_transcript_coords, exons_filt)

# map back to GRangesList, with group information
orf_in_genomic@unlistData$PID <- orf_ids[groupings(orf_in_genomic)]
#orf_in_genomic@unlistData$gene <- orf_gene_ids[groupings(orf_in_genomic)]

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
peptide_transcript_coords <- extract_peptide_coords(md, orf_transcript_coords_df)

# filter peptides for bad mappings
filtered_peptide_transcript_coords <- subset(peptide_transcript_coords, start(peptide_transcript_coords) != 0)
filtered_peptide_transcript_coords <- subset(filtered_peptide_transcript_coords, (mcols(filtered_peptide_transcript_coords)$txstart != mcols(filtered_peptide_transcript_coords)$pep_start))
filtered_peptide_transcript_coords <- subset(filtered_peptide_transcript_coords, (mcols(filtered_peptide_transcript_coords)$txend != mcols(filtered_peptide_transcript_coords)$pep_end))

# map to genomic coords
peptide_tx_names <- as.character(seqnames(filtered_peptide_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(filtered_peptide_transcript_coords) <- match(peptide_tx_names, names(exons_filt)) 

# original peptides in the GRanges
pep_ids <- filtered_peptide_transcript_coords$peptide
pep_PID_ids <- filtered_peptide_transcript_coords$PID
pep_gene_ids <- filtered_peptide_transcript_coords$gene_id

# ORFik map to genome coordinates
pep_in_genomic <- ORFik::pmapFromTranscriptF(filtered_peptide_transcript_coords, exons_filt, removeEmpty = F)

# map back to GRangesList, with group information
pep_in_genomic@unlistData$peptide <- pep_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$PID <- pep_PID_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$gene <- pep_gene_ids[groupings(pep_in_genomic)]

# remove 0 ranges
pep_in_genomic_gr <- unlist(pep_in_genomic, use.names=F) # convert to GRanges
# rename with transcript and peptide
tx_pep_names <- c(paste0(names(pep_in_genomic_gr), "_", pep_in_genomic_gr$peptide))
names(pep_in_genomic_gr) <- tx_pep_names # set names
# remove 0 ranges
pep_in_genomic_gr <- subset(pep_in_genomic_gr, (start(pep_in_genomic_gr) != 0 & end(pep_in_genomic_gr) != 0)  )

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(pep_in_genomic_gr), names(pep_in_genomic_gr), FUN = seq_along)
# add to GRanges
mcols(pep_in_genomic_gr)$exon_number <- exon_number_vec

# re-list
pep_in_genomic <- split(pep_in_genomic_gr, ~ names(pep_in_genomic_gr))

# ---------------------------------------- #


# ------------- export ------------- #

# PEPTIDES
# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, "peptides.bed12", rgb = 0)

# ORFS
# export bed12 of ORFs
ORFik::export.bed12(orf_in_genomic, "ORFs.bed12", rgb = 0)

# export GTF of ORFs
orf_in_genomic_gr$source <- c("custom")
orf_in_genomic_gr$type <- c("CDS")
orf_in_genomic_gr$phase <- 0
orf_in_genomic_gr$ORF_id <- names(orf_in_genomic_gr)
orf_in_genomic_gr$tx_id <- names(orf_in_genomic_gr)

export(orf_in_genomic_gr, "ORFs.gtf", format="gtf")

# TRANSCRIPTS
# export GTF of all transcripts that had mapped peptides
gtf_for_exporting <- import(gtf_import_file, format="gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]

export(gtf_filtered, "transcripts.gtf", format="gtf")

# ---------------------------------------- #


# ------- summary file -------- #

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
results_pept_df <- pep_in_genomic_gr %>% as_tibble()
results_pept_df <- separate(results_pept_df, txname, into = c("transcript_id", "peptide"), sep = "_", remove = TRUE)

# group by peptide and transcript to summarise based on how many exons peptide spans
results_pept_df_unique <- results_pept_df %>% 
  dplyr::group_by(peptide, transcript_id) %>% 
  dplyr::slice_max(exon_number) %>% 
  dplyr::mutate(number_exons = exon_number) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end, -width, -exon_number)

# merge results with metadata
metadata_to_merge <- md %>% 
  dplyr::select(PID, peptide, transcript_id, gene_id)

# merge results with metadata
# missing localisaiton etc
# metadata_to_merge <- md %>% 
#   dplyr::select(PID, peptide, transcript_id, gene_id, localisation, transcript_biotype, orf_type)

peptide_result <- merge(results_pept_df_unique, metadata_to_merge, by=c("PID", "peptide", "transcript_id"), all.x=T, all.y=F)

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

# define transcripts as known or novel
peptide_result <- peptide_result %>% 
  mutate(isoform_status = case_when(
    !startsWith(transcript_id, "EN") ~ "novel",
    TRUE ~ "known"
  ))

# determine ORF and peptide mapping status

# peptides mapped to a unique ORF (and gene) = high conf
# peptides mapped to multi ORFs from one gene = med conf
# peptides mapped to multi genomic locations = low conf

# NOTE some genes are the same with slightly different version numbers. Should we make it so script ignores versions??

peptide_result <- peptide_result %>% 
  dplyr::group_by(peptide) %>% 
  dplyr::mutate(pep_map_status = case_when(
    length(unique(gene_id))==1 & length(unique(PID))==1 ~ "high",
    length(unique(PID))>1 & length(unique(gene_id))==1 ~ "medium",
    length(unique(PID))>1 & length(unique(gene_id))>1 ~ "low",
    TRUE ~ NA)) %>% 
  ungroup()

# if an ORF is identified by at least one high confidence peptide, then the ORF is high confidence
peptide_result <- peptide_result %>% 
  dplyr::group_by(PID) %>% 
  dplyr::mutate(orf_map_status = case_when(
    pep_map_status %in% c("high") ~ "high",
    TRUE ~ "low")) %>% 
  dplyr::ungroup()

peptide_result <- peptide_result %>% 
  dplyr::group_by(peptide) %>% 
  dplyr::mutate(iso_map_status_1 = case_when(
    pep_map_status %in% c("high") & length(unique(transcript_id))==1 & !startsWith(transcript_id, "EN") ~ "novel_iso_uniq_identified",
    pep_map_status %in% c("high") & length(unique(transcript_id))==1 & startsWith(transcript_id, "EN") ~ "known_iso_uniq_identified",
    TRUE ~ NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(transcript_id) %>% 
  dplyr::mutate(iso_map_status = case_when(
    iso_map_status_1 %in% c("novel_iso_uniq_identified") ~ "novel_iso_uniq_identified",
    iso_map_status_1 %in% c("known_iso_uniq_identified") ~ "known_iso_uniq_identified",
    TRUE ~ "other")) %>% 
  dplyr::ungroup() %>% dplyr::select(-iso_map_status_1)

# summary(factor(peptide_result$single_mapped_protein))
# summary(factor(peptide_result$pep_map_status))
# summary(factor(peptide_result$orf_map_status))
# summary(factor(peptide_result$iso_map_status))

# include orf_map_status and pep_map_status in GTF mcols
results_pept_df$pep_map_status <- NULL
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
pep_in_genomic_gr_export$source <- c("custom")
pep_in_genomic_gr_export$type <- c("exon")

# export GTF of peptides
export(pep_in_genomic_gr_export, "peptides.gtf", format="gtf")

# export summary data
write.csv(peptide_result, "peptide_info.csv", row.names=F, quote=F)

# ---------------------------------------- #



