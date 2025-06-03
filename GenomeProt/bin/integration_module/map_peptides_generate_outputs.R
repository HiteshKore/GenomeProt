
# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

library(stringi)
source("global.R")
source("R/integration_functions.R")

option_list = list(
  make_option(c("-p", "--proteomics"), type="character", default=NULL,
              help="Proteomics data file", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Custom metadata used for proteomics", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF used to generate custom FASTA", metavar="character"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(scipen=999)

proteomics_import_file <- opt$proteomics
metadata_import_file <- opt$metadata
gtf_import_file <- opt$gtf
output_directory <- opt$savepath

# source("~/Documents/GenomeProt/GenomeProt/R/integration_functions.R")
# proteomics_import_file <- "/home/hiteshk/GenomeProt/PGData/Brain_data/peptide_Combined.txt"
# metadata_import_file <- "/home/hiteshk/GenomeProt/PGData/Brain_data/proteome_database_metadata.txt"
# gtf_import_file <- "/home/hiteshk/GenomeProt/PGData/Brain_data/proteome_database_transcripts.gtf"
# output_directory <- "/home/hiteshk/GenomeProt/PGData/Brain_data/"

# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))

gtf <- makeTxDbFromGFF(gtf_import_file, format="gtf") # make txdb of gtf

orf_df <- import_orf_metadata(metadata_import_file)

md<-integrate_metadata(pd,orf_df)


# ---------------------------------------- #


# ------------- run analysis ------------- #

# ------- map ORF transcript coords to spliced genomic coords ------- #

# create unique ID of ORF in transcript

md_order<-c("transcript_id","gene_id","PID","accession", "orf_genomic_coordinates", "chromosome","start","end","protein_sequence","protein_length","nt_length","tx_len","gene_name","strand","stranded_start","exon_start","exon_end","exon_rank","exon_length","total_length","txstart","txend")


md <- md %>%
  dplyr::select(any_of(md_order), everything())


md$orf_tx_id <- paste0(md$accession, "_", md$transcript_id)

# filter for unique orf and transcript for mapping the coordinates
orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id, strand)


orf_transcript_coords_df <- orf_transcript_coords_df[!(base::duplicated(orf_transcript_coords_df)),] # remove duplicates

# make GRanges from df of ORF transcript coordinates
orf_transcript_coords <- makeGRangesFromDataFrame(orf_transcript_coords_df,
                                                  keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                  seqnames.field="transcript_id", start.field="txstart", end.field="txend", strand.field="strand",
                                                  starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(orf_transcript_coords) <- c(orf_transcript_coords$orf_tx_id) # set names as unique IDs
mcols(orf_transcript_coords)$gene_id <- c(orf_transcript_coords_df$gene_id) # set gene_id

# get exons for mapping coordinates to genome
exons <- exonsBy(gtf, "tx", use.names=TRUE) # get exon data per transcript
exons_filt <- exons[names(exons) %in% orf_transcript_coords_df$transcript_id] # filter for only transcripts with mapped peptides

orf_tx_names <- as.character(seqnames(orf_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(orf_transcript_coords) <- match(orf_tx_names, names(exons_filt)) 

# create vector of unique names and gene ID for later
orf_ids <- orf_transcript_coords$orf_tx_id
orf_gene_ids <- orf_transcript_coords$gene_id

# use ORFik to map transcript coordinates to spliced genomic coordinates
orf_in_genomic <- ORFik::pmapFromTranscriptF(orf_transcript_coords, exons_filt, removeEmpty = T)

# map back to GRangesList, with group information, add gene_id back
orf_in_genomic@unlistData$PID <- orf_ids[groupings(orf_in_genomic)]
orf_in_genomic@unlistData$gene_id <- orf_gene_ids[groupings(orf_in_genomic)]

# unlist to add exon_number for GTF export
orf_in_genomic_gr <- unlist(orf_in_genomic, use.names=F) # convert to GRanges

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(orf_in_genomic_gr), mcols(orf_in_genomic_gr)$PID, FUN = seq_along)

# add to GRanges
mcols(orf_in_genomic_gr)$exon_number <- exon_number_vec
# re-list
orf_in_genomic <- split(orf_in_genomic_gr, ~ mcols(orf_in_genomic_gr)$PID)

# ---------------------------------------- #


# ------- map peptide transcript coords to spliced genomic coords ------- #

# use ORF transcript coords to determine peptide transcript coords
peptide_transcript_coords <- extract_peptide_coords(md)

# get vector of transcript names
peptide_tx_names <- as.character(seqnames(peptide_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(peptide_transcript_coords) <- match(peptide_tx_names, names(exons_filt)) 

# create vectors of unique names and gene ID for later
pep_ids <- peptide_transcript_coords$peptide
pep_PID_ids <- peptide_transcript_coords$PID
pep_gene_ids <- peptide_transcript_coords$gene_id

# use ORFik to map transcript coordinates to spliced genomic coordinates
pep_in_genomic <- ORFik::pmapFromTranscriptF(peptide_transcript_coords, exons_filt, removeEmpty = F)

# map back to GRangesList, with group information, add other info back
pep_in_genomic@unlistData$peptide <- pep_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$PID <- pep_PID_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$gene <- pep_gene_ids[groupings(pep_in_genomic)]

# unlist to add exon_number for GTF export
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


# ------------- export BED12 files ------------- #

# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, paste0(output_directory, "/peptides.bed12"), rgb = 0)

# export bed12 of ORFs
# currently export with PID_transcript as the name, means ORFs are often present multiple times
ORFik::export.bed12(orf_in_genomic, paste0(output_directory, "/ORFs.bed12"), rgb = 0)

# format GTF of all transcripts that had mapped peptides
gtf_for_exporting <- rtracklayer::import(gtf_import_file, format="gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript_id]
gtf_filtered$group_id <- "transcripts"

# reformat exons for bed12
gtf_as_bed12 <- gtf_filtered[mcols(gtf_filtered)$type == "exon"]

names(gtf_as_bed12) <- paste0(gtf_as_bed12$transcript_id, "_", gtf_as_bed12$gene_id)

# convert to grl
tx_in_genomic <- split(gtf_as_bed12, ~ names(gtf_as_bed12))

# export bed12 of transcripts
ORFik::export.bed12(tx_in_genomic, paste0(output_directory, "/transcripts.bed12"), rgb = 0)

# ---------------------------------------- #

# ------- transcripts and IsoVis GTF -------- #

# for combined GTF

# create object unformatted for IsoVis
orf_in_genomic_gr_isovis <- orf_in_genomic_gr

# format GTF to combine with transcripts and peptides last
orf_in_genomic_gr$source <- c("custom")
orf_in_genomic_gr$type <- c("CDS")
orf_in_genomic_gr$phase <- 0
orf_in_genomic_gr$ORF_id <- names(orf_in_genomic_gr)
orf_in_genomic_gr$transcript_id <- names(orf_in_genomic_gr)

names(orf_in_genomic_gr) <- NULL
orf_in_genomic_gr$group_id <- "ORFs"

# for IsoVis

# re-filter initial GTF of transcripts and exons
gtf_isovis <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]

# set standard GTF columns
orf_in_genomic_gr_isovis$source <- c("custom")
orf_in_genomic_gr_isovis$type <- c("CDS")
orf_in_genomic_gr_isovis$phase <- 0
orf_in_genomic_gr_isovis$transcript_id <- names(orf_in_genomic_gr)

# get PID_transcript column as df
pids <- data.frame(ids = orf_in_genomic_gr_isovis$PID)
# separate into just protein ID based on last occurrence of an '_'
pids <- pids %>% separate(ids, into = "protein_id", sep = "\\_(?!.*_)", remove=F)
# paste as protein_id
orf_in_genomic_gr_isovis$protein_id <- paste0(pids$protein_id)
# remove other columns
names(orf_in_genomic_gr_isovis) <- NULL
orf_in_genomic_gr_isovis$PID <- NULL

# combine transcripts, exons and CDS
isovis_export <- c(gtf_isovis, orf_in_genomic_gr_isovis)
# sort
isovis_export <- sortSeqlevels(isovis_export)
isovis_export <- sort(isovis_export)
# export GTF compatible with IsoVis
rtracklayer::export(isovis_export, paste0(output_directory, "/transcripts_and_ORFs_for_isovis.gtf"), format="gtf")

# ---------------------------------------- #

# ------- summary file of peptide mappings -------- #

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
results_pept_df <- pep_in_genomic_gr %>% as_tibble()
results_pept_df <- separate(results_pept_df, txname, into = c("transcript_id", "peptide"), sep = "_", remove = TRUE)
results_pept_df$gene_id <- results_pept_df$gene
results_pept_df$gene <- NULL

# group by peptide and transcript to summarise based on how many exons peptide spans
results_pept_df_unique <- results_pept_df %>% 
  dplyr::group_by(PID, peptide, transcript_id) %>% 
  dplyr::slice_max(exon_number) %>% 
  dplyr::mutate(number_exons = exon_number) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end, -width, -exon_number,-strand)

# merge results with metadata
 #metadata_to_merge <- md %>% 
 #  dplyr::select(PID, peptide, transcript_id, gene_id, gene_name, protein_length, tx_len)
 
 peptide_result <- left_join(results_pept_df_unique, md, by=c("PID", "peptide", "gene_id", "transcript_id"))
 
 peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]
 
 
 peptide_result$seqnames <- NULL
 
 peptide_result <- peptide_result %>% 
   dplyr::mutate(longest_orf_in_transcript = case_when(
     longest_orf_in_transcript == "Y" ~ TRUE,
     longest_orf_in_transcript == "N" ~ FALSE
   ))

setDT(peptide_result)

peptide_result[, c("peptide_ids_gene", "peptide_ids_orf", "peptide_ids_transcript", "shared_novel_protein_peptide") := 
                 .(length(unique(gene_id)) == 1,
                   length(unique(PID)) == 1,
                   length(unique(transcript_id)) == 1,
                   length(unique(PID)) > 1 & all(startsWith(PID, "ORF"))),
               by = peptide]

peptide_result[, orf_identified := any(peptide_ids_orf == TRUE), by = PID]
peptide_result[, gene_identified := any(peptide_ids_gene == TRUE), by = gene_id]
peptide_result[, transcript_identified := any(peptide_ids_transcript == TRUE), by = transcript_id]

# get missing peptides back in the output
missing_peptides <- pd %>% dplyr::filter(!(peptide %in% peptide_result$peptide))

# ensure same col names
columns_to_add <- setdiff(names(peptide_result), names(missing_peptides))

# add missing columns to missing df with NA values
missing_peptides[columns_to_add] <- NA

# ensure the column order is same before rbind
missing_peptides <- missing_peptides[, names(peptide_result)]

combined_peptide_result <- rbind(peptide_result, missing_peptides)

combined_peptide_result$PID <- gsub(",", ".", combined_peptide_result$PID)
combined_peptide_result <- combined_peptide_result[!(base::duplicated(combined_peptide_result)),]
combined_peptide_result$transcript_length <- combined_peptide_result$tx_len

combined_peptide_result <- combined_peptide_result %>% 
  mutate(simplified_biotype = case_when(
    transcript_biotype %in% c("protein_coding", "protein_coding_LoF", "protein_coding_CDS_not_defined") ~ "protein_coding",
    transcript_biotype %in% c("polymorphic_pseudogene", "pseudogene", "processed_pseudogene", 
                              "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene",
                              "unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_processed_pseudogene",
                              "translated_unprocessed_pseudogene", "unitary_pseudogene") ~ "pseudogene",
    transcript_biotype %in% c("nonsense_mediated_decay", "non_stop_decay") ~ "NMD",
    transcript_biotype %in% c("retained_intron") ~ "retained_intron",
    transcript_biotype %in% c("lncRNA") ~ "lncRNA",
    transcript_biotype %in% c("novel") ~ "novel",
    TRUE ~ "other"
  ))

# rearrange columns for output
#combined_peptide_result<-combined_peptide_result%>%rename(orf_genomic_coordinates=location)

combined_peptide_result <- combined_peptide_result %>%
  group_by(accession) %>%
  mutate(is_unique_location = n_distinct(orf_genomic_coordinates) == 1) %>%
  ungroup()


combined_peptide_result <- combined_peptide_result %>% dplyr::select(peptide,accession,PID,transcript_id,gene_id,gene_name,strand,
                                                                     number_exons,transcript_length,transcript_biotype,simplified_biotype,
                                                                     protein_length,orf_genomic_coordinates,orf_type,localisation,uniprot_status,
                                                                     openprot_id,`molecular_weight(kDA)`,isoelectric_point,hydrophobicity,
                                                                     aliphatic_index,longest_orf_in_transcript,is_unique_location,peptide_ids_gene,peptide_ids_orf,
                                                                     peptide_ids_transcript,shared_novel_protein_peptide,orf_identified,
                                                                     gene_identified,transcript_identified)

# export summary data
write.csv(combined_peptide_result, paste0(output_directory, "/peptide_info.csv"), row.names=F, quote=F)

# ---------------------------------------- #


# ------- combined GTF -------- #

# include orf_status and peptide_status in GTF mcols

peptide_result<-peptide_result%>%dplyr::select("PID","transcript_id","peptide","gene_id","number_exons","strand" ,"gene_name","protein_length","tx_len","accession","orf_genomic_coordinates","transcript_biotype","orf_type","localisation","uniprot_status",
"longest_orf_in_transcript","peptide_ids_gene","peptide_ids_orf" ,"peptide_ids_transcript", "shared_novel_protein_peptide","orf_identified","gene_identified","transcript_identified")
                        

results_to_merge_with_granges <- left_join(results_pept_df, peptide_result, by=c("transcript_id", "peptide", "strand", "PID", "gene_id"))
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

# export annotations for vis
combined <- c(pep_in_genomic_gr_export, orf_in_genomic_gr, gtf_filtered)
rtracklayer::export(combined, paste0(output_directory, "/combined_annotations.gtf"), format="gtf")

# ---------------------------------------- #



