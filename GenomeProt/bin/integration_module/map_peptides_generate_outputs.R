# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

# functions

# import proteomics file (fragpipe, maxquant etc)
import_proteomics_data <- function(proteomics_file) {
  # required columns
  # Peptide, Protein, Mapped Proteins
  # optional: Protein Start
  
  # import proteomics
  protfile <- data.table::fread(proteomics_file)

  # ensure compatibility
  if (all(c("Peptide", "Protein", "Mapped_proteins") %in% colnames(protfile))) {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, Mapped_proteins)
    protfile <- protfile %>% dplyr::mutate(
      all_mappings = dplyr::case_when(
        (is.na(Mapped_proteins) == TRUE | nchar(Mapped_proteins) == 0) ~ paste0(Protein),
        TRUE ~ paste0(Protein, ", ", Mapped_proteins)
      )
    )
  }

  # separate into one row per mapped ORF
  prot_expanded <- tidyr::separate_rows(protfile, all_mappings, sep = ",")

  # remove white space and characters
  prot_expanded <- prot_expanded %>% 
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ gsub(" ", "", .))) %>%
    dplyr::mutate(accession = all_mappings)

  # remove any duplications
  prot_expanded <- prot_expanded[!(base::duplicated(prot_expanded)),]

  # select columns to keep
  prot_expanded <- prot_expanded %>% dplyr::select(accession, Peptide)

  # remove peptides mapped to reversed sequences
  prot_expanded <- prot_expanded %>% 
    dplyr::filter(!startsWith(accession, "rev"))

  # rename column
  prot_expanded$peptide <- prot_expanded$Peptide
  prot_expanded$Peptide <- NULL

  return(prot_expanded)  
}

###
# import custom database metadata of ORFs
import_orf_metadata <- function(metadata_file) {
  orf_data <- data.table::fread(metadata_file)
  orf_data$transcript_id <- orf_data$transcript
  orf_data$transcript <- NULL
  orf_data$PID <- paste0(orf_data$accession, "|CO=", orf_data$orf_genomic_coordinates)
  orf_data <- orf_data %>% dplyr::select(PID, accession, protein_description, orf_genomic_coordinates, gene, gene_symbol, transcript_id, strand, transcript_biotype, transcript_coordinates, orf_type, localisation, uniprot_status, openprot_id,
                                         protein_sequence, `molecular_weight(kDA)`, isoelectric_point, hydrophobicity, aliphatic_index, longest_orf_in_transcript)
  return(orf_data)
}

# import custom database FASTA of ORFs


integrate_metadata<-function(pd,orf_df){
  #subset metadata for proteins detected in proteomics

  orf_df_protein_detected <- orf_df %>% dplyr::filter(accession %in% pd$accession) %>%
    dplyr::rename(gene_id = gene, gene_name = gene_symbol) %>%
    tidyr::separate(orf_genomic_coordinates, into = c("chromosome", "start", "end"), sep = "\\:|-", remove = FALSE) %>%
    dplyr::mutate(protein_length = nchar(protein_sequence))

  # get transcript lengths
  gtf_txdb <- txdbmaker::makeTxDbFromGFF(gtf_import_file, format = "gtf")
  tx_lengths <- transcriptLengths(gtf_txdb)
  tx_lengths$transcript_id <- tx_lengths$tx_name
  tx_lengths <- tx_lengths %>% dplyr::select(transcript_id, tx_len)
  
  # get gene names
  gtf_import <- rtracklayer::import(gtf_import_file, format = "gtf") %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(type == "transcript") %>% 
    dplyr::select(transcript_id, gene_name, strand)

  # merge tx info
  tx_data <- merge(tx_lengths, gtf_import, by = "transcript_id", all.x = T, all.y = F) # transcript_id tx_len gene_name strand
  # get transcripts for mapping
  txs <- exonsBy(gtf_txdb, by = c("tx"), use.names = T)

  orf_df_protein_detected <- dplyr::left_join(orf_df_protein_detected, tx_lengths, by = "transcript_id")

  orf_df_protein_detected$start <- as.numeric(orf_df_protein_detected$start)
  orf_df_protein_detected$end <- as.numeric(orf_df_protein_detected$end)

  # re orient start and end to match strand for mapping
  orf_df_protein_detected <- orf_df_protein_detected %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(stranded_start = dplyr::case_when(
      strand == "+" ~ min(start,end),
      strand == "-" ~ max(start,end),
    ))

  txs_filtered <- txs[names(txs) %in% orf_df_protein_detected$transcript_id]
  txs_unlisted <- unlist(txs_filtered)
  # function to convert genomic start position to transcript position
  convert_gene_pos_to_transcript_pos <- function(txdb, input_df) {
    # get all relevant exons
    all_exons <- txdb[names(txdb) %in% input_df$transcript_id]

    exon_data <- data.frame(
      transcript_id = names(all_exons),
      exon_start = start(all_exons),
      exon_end = end(all_exons),
      exon_rank = mcols(all_exons)$exon_rank
    ) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(exon_rank) %>%
      dplyr::mutate(
        exon_length = exon_end - exon_start + 1,
        total_length = cumsum(exon_length)
      )

    # join input data with exon data
    result <- input_df %>%
      dplyr::left_join(exon_data, by = "transcript_id", relationship = "many-to-many") %>%
      dplyr::filter(stranded_start >= exon_start & stranded_start <= exon_end) %>%
      dplyr::group_by(accession,transcript_id, stranded_start) %>%
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(
        txstart = dplyr::case_when(
          strand == "+" & exon_rank == 1 ~ stranded_start - exon_start,
          strand == "+" ~ (stranded_start - exon_start) + (total_length - exon_length),
          strand == "-" & exon_rank == 1 ~ exon_end - stranded_start,
          strand == "-" ~ (exon_end - stranded_start) + (total_length - exon_length)
        )
      )
    return(result)
  }

  # apply to orf genomic coords
  orf_transcriptomic_coords <- convert_gene_pos_to_transcript_pos(txs_unlisted, orf_df_protein_detected)

  # ORFs that could not be mapped to transcripts are returned with -1 starts
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txstart >= 0)
  orf_transcriptomic_coords$nt_length <- orf_transcriptomic_coords$protein_length * 3
  # get transcript end coord
  orf_transcriptomic_coords$txend <- orf_transcriptomic_coords$txstart + orf_transcriptomic_coords$nt_length - 1

  # remove ORF coords outside transcript ends
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txend < tx_len)

  # first filter proteomics data to remove peptides mapping to ORFs with multiple loci
  #multi_loci_peptides <- proteomics_data[grepl("\\,chr", proteomics_data$PID),]

  #proteomics_filtered <- proteomics_data %>% dplyr::filter(!(peptide %in% multi_loci_peptides$peptide))

  # combine with proteomics data
  proteome_data <- pd %>% dplyr::select(accession, peptide) %>% base::unique()

  metadata <- dplyr::full_join(orf_transcriptomic_coords, proteome_data, by = "accession")
  metadata <- metadata[!(base::duplicated(metadata)),]
  metadata <- metadata %>% dplyr::filter(!is.na(transcript_id) & !is.na(txstart))

  find_peptide_position <- function(peptides, protein_sequences) {
    # exact match
    exact_matches <- stringi::stri_locate_first_fixed(protein_sequences, peptides)

    # non-exact matches use fuzzy matching
    non_exact_mask <- is.na(exact_matches[,1])

    if (any(non_exact_mask)) {
      fuzzy_matches <- vapply(which(non_exact_mask), function(i) {
        peptide <- peptides[i]
        protein <- protein_sequences[i]

        # generate regex pattern allowing up to 3 mismatches
        pattern <- paste0(strsplit(peptide, "")[[1]], collapse = "(.{0,1})")

        match <- stringi::stri_locate_first_regex(protein, pattern)
        if (is.na(match[1])) -1L else as.integer(match[1])
      }, integer(1))

      exact_matches[non_exact_mask,1] <- fuzzy_matches
    }
    return(as.integer(exact_matches[,1]))
  }

  # apply to dataframe
  data.table::setDT(metadata)
  metadata[, mapped_pep_start := find_peptide_position(peptide, protein_sequence)]

  # if exact match or one mismatch isn't found, use proteomics defined start site
  if (c("Protein Start") %in% colnames(metadata)) {
    metadata <- metadata %>% dplyr::mutate(pep_start = dplyr::case_when(
      (mapped_pep_start != `Protein Start` & mapped_pep_start != -1) ~ mapped_pep_start,
      (mapped_pep_start != `Protein Start` & mapped_pep_start == -1) ~ `Protein Start`,
      TRUE ~ mapped_pep_start
    ))
  } else {
    metadata <- metadata %>% dplyr::filter(mapped_pep_start != -1 & !is.na(mapped_pep_start))
    metadata$pep_start <- metadata$mapped_pep_start
  }

  # get peptide end location within every ORF
  metadata$pep_end <- metadata$pep_start + nchar(metadata$peptide) - 1

  metadata_output <- metadata %>% dplyr::filter(pep_start < protein_length & pep_end <= protein_length)
  metadata_output$mapped_pep_start <- NULL

  return(metadata)
}

# get peptide coords from transcript coords
extract_peptide_coords <- function(metadata_df) {
  # subset metadata for conversion of peptide coords
  txcoordsdf_subset <- metadata_df %>% dplyr::select(orf_tx_id, transcript_id, PID, gene_id, pep_start, pep_end, peptide, txstart, txend)

  # get peptide locations within every transcript
  txcoordsdf_subset$txstart <- as.numeric(txcoordsdf_subset$txstart)
  txcoordsdf_subset$txend <- as.numeric(txcoordsdf_subset$txend)

  txcoordsdf_subset$start <- txcoordsdf_subset$txstart + ((txcoordsdf_subset$pep_start-1) * 3)
  txcoordsdf_subset$end <- txcoordsdf_subset$start + (nchar(txcoordsdf_subset$peptide) * 3)
  txcoordsdf_subset$peptide <- as.factor(txcoordsdf_subset$peptide)

  # sometimes peptide ends are off by 1, if so correct them to transcript end site
  txcoordsdf_subset <- txcoordsdf_subset %>% 
    dplyr::mutate(end_corrected = dplyr::case_when(
      (end == txend + 1) ~ txend,
      TRUE ~ end
    ))

  # replace original end col
  txcoordsdf_subset$end <- txcoordsdf_subset$end_corrected

  # remove tmp end col
  txcoordsdf_subset$end_corrected <- NULL

  # ensure peptides are within transcript coords
  txcoordsdf_subset <- txcoordsdf_subset %>% 
    dplyr::filter(start >= txstart & end <= txend)

  # GRanges from df of transcriptomic coords of peptides
  coords_peptides <- GenomicRanges::makeGRangesFromDataFrame(txcoordsdf_subset,
                                                             keep.extra.columns = TRUE,
                                                             ignore.strand = FALSE,
                                                             seqinfo = NULL,
                                                             seqnames.field = "transcript_id",
                                                             start.field = "start",
                                                             end.field = "end",
                                                             strand.field = "strand",
                                                             starts.in.df.are.0based = TRUE,
                                                             na.rm = TRUE)
  # set names
  names(coords_peptides) <- c(coords_peptides$transcript)

  # filter peptides for bad mappings
  coords_peptides <- subset(coords_peptides, start(coords_peptides) != 0)
  coords_peptides <- subset(coords_peptides, end(coords_peptides) <= coords_peptides$txend)

  return(coords_peptides)
}

suppressPackageStartupMessages({
  library(optparse)
})

# define options
option_list = list(
  make_option(c("-p", "--proteomics"), type="character", default=NULL,
              help="Reformatted proteomics results file (e.g. peptide_data.tsv)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Proteome database metadata file (e.g. proteome_database_metadata.txt)", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Proteome database transcripts GTF file (e.g. proteome_database_transcripts.gtf)", metavar="character"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# store the input arguments
proteomics_import_file <- opt$proteomics
metadata_import_file <- opt$metadata
gtf_import_file <- opt$gtf
output_directory <- opt$savepath

# check if any input arguments are missing
if (is.null(proteomics_import_file)) {
  stop("Please provide a reformatted proteomics results file.")
} else if (is.null(metadata_import_file)) {
  stop("Please provide a proteome database metadata file.")
} else if (is.null(gtf_import_file)) {
  stop("Please provide a proteome database transcripts GTF file.")
} else if (is.null(output_directory)) {
  stop("Please specify an output directory.")
}

# check the input arguments
if (!file.exists(proteomics_import_file)) {
  stop(paste0("The reformatted proteomics results file '", proteomics_import_file, "' does not exist."))
} else if (!file.exists(metadata_import_file)) {
  stop(paste0("The proteome database metadata file '", metadata_import_file, "' does not exist."))
} else if (!file.exists(gtf_import_file)) {
  stop(paste0("The proteome database transcripts GTF file '", gtf_import_file, "' does not exist."))
} else if (file.exists(output_directory) && !dir.exists(output_directory)) {
  stop(paste0("'", output_directory, "' exists but is not a directory."))
}

# ensure that none of the specified files are empty
if (file.size(proteomics_import_file) == 0) {
  stop(paste0("The reformatted proteomics results file '", proteomics_import_file, "' must not be empty."))
} else if (file.size(metadata_import_file) == 0) {
  stop(paste0("The proteome database metadata file '", metadata_import_file, "' must not be empty."))
} else if (file.size(gtf_import_file) == 0) {
  stop(paste0("The proteome database transcripts GTF file '", gtf_import_file, "' must not be empty."))
}

# create the output directory if it does not already exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

options(scipen = 999)

# load the rest of the libraries needed for the remainder of the script to run
suppressPackageStartupMessages({
  library(txdbmaker)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(ORFik)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(stringi)
  library(stats)
})

# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))
gtf <- txdbmaker::makeTxDbFromGFF(gtf_import_file, format = "gtf")
orf_df <- import_orf_metadata(metadata_import_file)
md <- integrate_metadata(pd, orf_df)

# ------------- run analysis ------------- #

# ------- map ORF transcript coords to spliced genomic coords ------- #

# create unique ID of ORF in transcript

md_order <- c("transcript_id", "gene_id", "PID", "accession", "orf_genomic_coordinates", "chromosome", "start", "end", "protein_sequence", "protein_length", "nt_length", "tx_len", "gene_name", "strand", "stranded_start", "exon_start", "exon_end", "exon_rank", "exon_length", "total_length", "txstart", "txend")

md <- md %>%
  dplyr::select(dplyr::any_of(md_order), dplyr::everything())

md$orf_tx_id <- paste0(md$accession, "_", md$transcript_id)

# filter for unique orf and transcript for mapping the coordinates
orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id, strand)

orf_transcript_coords_df <- orf_transcript_coords_df[!(base::duplicated(orf_transcript_coords_df)),] # remove duplicates

# make GRanges from df of ORF transcript coordinates
orf_transcript_coords <- GenomicRanges::makeGRangesFromDataFrame(orf_transcript_coords_df,
                                                                 keep.extra.columns = TRUE, ignore.strand = FALSE, seqinfo = NULL,
                                                                 seqnames.field = "transcript_id", start.field = "txstart", end.field = "txend", strand.field = "strand",
                                                                 starts.in.df.are.0based = FALSE, na.rm = TRUE)

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
exon_number_vec <- stats::ave(seq_along(orf_in_genomic_gr), mcols(orf_in_genomic_gr)$PID, FUN = seq_along)

# add to GRanges
mcols(orf_in_genomic_gr)$exon_number <- exon_number_vec
# re-list
orf_in_genomic <- split(orf_in_genomic_gr, ~ mcols(orf_in_genomic_gr)$PID)

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
pep_in_genomic_gr$tx_PID_grouping <- paste0(pep_in_genomic_gr$PID, "_", names(pep_in_genomic_gr))

# remove 0 ranges
pep_in_genomic_gr <- subset(pep_in_genomic_gr, (start(pep_in_genomic_gr) != 0 & end(pep_in_genomic_gr) != 0)  )

# create vector of exon number per peptide and transcript
exon_number_vec <- stats::ave(seq_along(pep_in_genomic_gr), pep_in_genomic_gr$tx_PID_grouping, FUN = seq_along)

# add to GRanges
mcols(pep_in_genomic_gr)$exon_number <- exon_number_vec

mcols(pep_in_genomic_gr)$tx_pid_grouping <- NULL

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
mcols(pep_in_genomic_gr)$transcript_id <- names(pep_in_genomic_gr)
mcols(pep_in_genomic_gr)$transcript_id <- word(names(pep_in_genomic_gr),1,sep="_")
mcols(pep_in_genomic_gr)$peptide <- word(names(pep_in_genomic_gr),2,sep="_")

results_pept_df <- pep_in_genomic_gr %>% tibble::as_tibble()
results_pept_df$gene_id <- results_pept_df$gene
results_pept_df$gene <- NULL

# extract transcript exon overlap with peptide coordinates
exons_gr<-unlist(exons_filt,use.names = TRUE)
mcols(exons_gr)$transcript_id<-names(exons_gr)

hits<-findOverlaps(pep_in_genomic_gr,exons_gr)

overlapping_peptides<-pep_in_genomic_gr[queryHits(hits)]
overlapping_exons<-exons_gr[subjectHits(hits)]

# Combine information into a DataFrame
peptide_exon_overlap <- data.frame(
  peptide_coordinates_in_exon = paste0(seqnames(overlapping_peptides),":",start(overlapping_peptides),"-",end(overlapping_peptides)),
  strand = strand(overlapping_peptides),
  transcript_id = overlapping_peptides$transcript_id,
  gene_id = overlapping_peptides$gene,
  exon_number = overlapping_peptides$exon_number,
  peptide = overlapping_peptides$peptide,
  PID = overlapping_peptides$PID,
  peptide_overlapping_exon_coordinates = paste0(seqnames(overlapping_exons),":",start(overlapping_exons),"-",end(overlapping_exons)),
  exon_transcript = overlapping_exons$transcript_id,
  exon_number_transcript = overlapping_exons$exon_rank
) %>% dplyr::mutate(
  tx_PID_grouping = paste0(PID,"_",transcript_id,"_",peptide),
)

overlap_df_reduced <- peptide_exon_overlap %>% dplyr::filter(transcript_id==exon_transcript) %>% base::unique() %>%
  dplyr::group_by(tx_PID_grouping) %>%
  dplyr::summarise(
    overlapping_peptide_coordinates = paste0(base::unique(peptide_coordinates_in_exon),collapse =","),
    overlapping_exon_coordinates = paste0(base::unique(peptide_overlapping_exon_coordinates),collapse =","),
    overlapping_exon_number = paste0(base::unique(exon_number_transcript),collapse =","),
    number_of_mapped_exons = dplyr::n_distinct(peptide_overlapping_exon_coordinates),
    strand_transcript = dplyr::first(strand),
    peptide = dplyr::first(peptide),
    PID = dplyr::first(PID),
    gene_id = dplyr::first(gene_id),
    transcript_id = dplyr::first(transcript_id),
    mapped_exons = paste(base::unique(exon_number),collapse=","),
    .groups = 'drop'
  )

overlap_df_reduced <- overlap_df_reduced %>%
  dplyr::mutate(
    peptide_genomic_coordinates = purrr::map2_chr(overlapping_peptide_coordinates, strand_transcript, function(coord_str, strand_info) {
      coords <- stringr::str_split(coord_str, ",")[[1]]
      if (strand_info == "+") {
        first_coord <- coords[1]
        last_coord <- coords[length(coords)]
      } else {
        first_coord <- coords[length(coords)]
        last_coord <- coords[1]
      }
      chr_first <- stringr::str_extract(first_coord, "^[^:]+")
      start_first <- stringr::str_extract(first_coord, "(?<=:)[0-9]+")
      end_last <- stringr::str_extract(last_coord, "(?<=-)[0-9]+")
      paste0(chr_first, ":", start_first, "-", end_last)
    })
  ) %>%
  dplyr::mutate(
    peptide_mapping_type = dplyr::case_when(
      lengths(strsplit(mapped_exons,",")) >1 ~ "exon-spanning peptide",
      TRUE~ "exonic peptide"
    )
  )

mcols(pep_in_genomic_gr)$tx_PID_grouping <- NULL

# re-list
pep_in_genomic <- split(pep_in_genomic_gr, ~ names(pep_in_genomic_gr))

# ------------- export BED12 files ------------- #

# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, file.path(output_directory, "peptides.bed12"), rgb = 0)

# export bed12 of ORFs
# currently export with PID_transcript as the name, means ORFs are often present multiple times
ORFik::export.bed12(orf_in_genomic, file.path(output_directory, "ORFs.bed12"), rgb = 0)

# format GTF of all transcripts that had mapped peptides
gtf_for_exporting <- rtracklayer::import(gtf_import_file, format = "gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript_id]
gtf_filtered$group_id <- "transcripts"

# reformat exons for bed12
gtf_as_bed12 <- gtf_filtered[mcols(gtf_filtered)$type == "exon"]

names(gtf_as_bed12) <- paste0(gtf_as_bed12$transcript_id, "_", gtf_as_bed12$gene_id)

# convert to grl
tx_in_genomic <- split(gtf_as_bed12, ~ names(gtf_as_bed12))

# export bed12 of transcripts
ORFik::export.bed12(tx_in_genomic, file.path(output_directory, "transcripts.bed12"), rgb = 0)

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
PIDs <- data.frame(ids = orf_in_genomic_gr_isovis$PID)
# separate into just protein ID based on last occurrence of an '_'
PIDs <- PIDs %>% tidyr::separate(ids, into = "protein_id", sep = "\\_(?!.*_)", remove=F)
# paste as protein_id
orf_in_genomic_gr_isovis$protein_id <- paste0(PIDs$protein_id)
# remove other columns
names(orf_in_genomic_gr_isovis) <- NULL
orf_in_genomic_gr_isovis$PID <- NULL

# combine transcripts, exons and CDS
isovis_export <- c(gtf_isovis, orf_in_genomic_gr_isovis)
# sort
isovis_export <- Seqinfo::sortSeqlevels(isovis_export)
isovis_export <- sort(isovis_export)
# export GTF compatible with IsoVis
rtracklayer::export(isovis_export, file.path(output_directory, "transcripts_and_ORFs_for_isovis.gtf"), format = "gtf")

# ------- summary file of peptide mappings -------- #

# merge results with metadata

peptide_result <- dplyr::left_join(overlap_df_reduced, md, by=c("PID", "peptide", "gene_id", "transcript_id"))

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

peptide_result <- peptide_result %>%
  dplyr::mutate(longest_orf_in_transcript = dplyr::case_when(
    longest_orf_in_transcript == "Y" ~ TRUE,
    longest_orf_in_transcript == "N" ~ FALSE
  ))

data.table::setDT(peptide_result)

peptide_result[, c("peptide_ids_gene", "peptide_ids_orf", "peptide_ids_transcript", "shared_novel_protein_peptide") :=
                 .(length(base::unique(gene_id)) == 1,
                   length(base::unique(accession)) == 1,
                   length(base::unique(transcript_id)) == 1,
                   length(base::unique(accession)) > 1 & all(startsWith(accession, "ORF"))),
               by = peptide]

peptide_result[, orf_identified := any(peptide_ids_orf == TRUE), by = accession]
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
  dplyr::mutate(simplified_biotype = dplyr::case_when(
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
combined_peptide_result <- combined_peptide_result %>%
  dplyr::group_by(accession) %>%
  dplyr::mutate(is_unique_location = dplyr::n_distinct(orf_genomic_coordinates) == 1) %>%
  dplyr::ungroup()

combined_peptide_result <- combined_peptide_result %>% dplyr::select(peptide,peptide_genomic_coordinates,accession,PID,protein_description,transcript_id,gene_id,gene_name,strand,
                                                                     overlapping_peptide_coordinates,overlapping_exon_coordinates,overlapping_exon_number,number_of_mapped_exons,peptide_mapping_type,transcript_length,transcript_biotype,simplified_biotype,
                                                                     protein_length,orf_genomic_coordinates,orf_type,localisation,uniprot_status,
                                                                     openprot_id,`molecular_weight(kDA)`,isoelectric_point,hydrophobicity,
                                                                     aliphatic_index,longest_orf_in_transcript,is_unique_location,peptide_ids_gene,peptide_ids_orf,
                                                                     peptide_ids_transcript,shared_novel_protein_peptide,orf_identified,
                                                                     gene_identified,transcript_identified)

# export summary data
readr::write_tsv(combined_peptide_result, file.path(output_directory, "peptide_info.tsv"))

# ------- combined GTF -------- #

# include orf_status and peptide_status in GTF mcols

peptide_result<-peptide_result%>%dplyr::select("PID","transcript_id","peptide","gene_id","number_of_mapped_exons","strand" ,"gene_name","protein_length","tx_len","accession","orf_genomic_coordinates","transcript_biotype","orf_type","localisation","uniprot_status",
"longest_orf_in_transcript","peptide_ids_gene","peptide_ids_orf" ,"peptide_ids_transcript", "shared_novel_protein_peptide","orf_identified","gene_identified","transcript_identified")

results_to_merge_with_granges <- dplyr::left_join(results_pept_df, peptide_result, by=c("transcript_id", "peptide", "strand", "PID", "gene_id"))
results_to_merge_with_granges <- results_to_merge_with_granges[!(duplicated(results_to_merge_with_granges)),]

results_to_merge_with_granges$naming <- paste0(results_to_merge_with_granges$transcript_id, "_", results_to_merge_with_granges$peptide)

# make GRanges from df of ORF transcript coordinates
pep_in_genomic_gr_export <- GenomicRanges::makeGRangesFromDataFrame(results_to_merge_with_granges,
                                                                    keep.extra.columns = TRUE, ignore.strand = FALSE, seqinfo = NULL,
                                                                    seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand",
                                                                    starts.in.df.are.0based = FALSE, na.rm = TRUE)

names(pep_in_genomic_gr_export) <- c(pep_in_genomic_gr_export$naming) # set names

# add mcols
pep_in_genomic_gr_export$source <- "custom"
pep_in_genomic_gr_export$type <- "exon"
pep_in_genomic_gr_export$group_id <- "peptides"

# export annotations for vis
combined <- c(pep_in_genomic_gr_export, orf_in_genomic_gr, gtf_filtered)

gtf_file <- file.path(output_directory, "combined_annotations.gtf")
rtracklayer::export(combined, gtf_file, format = "gtf")

# read exported file
gtf_lines <- readLines(gtf_file)
custom_header <- c("##GenomeProt", "##transcript,ORF and peptide annotations")
writeLines(c(custom_header, gtf_lines), gtf_file)