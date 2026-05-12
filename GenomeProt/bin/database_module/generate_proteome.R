# import gtf and filter for minimum transcript counts
filter_custom_gtf <- function(customgtf, organism, tx_counts = NA, min_count = NA, outdir) {
  # import bambu gtf
  bambu_data <- rtracklayer::import(customgtf)

  if (!missing(tx_counts)) { # if counts file has been uploaded
    # define default value
    if (missing(min_count)) {
      min_count <- 5
    }

    # read in counts
    counts <- fread(tx_counts)

    # filter for counts greater than or equal to min_count in any numeric column
    counts_filt <- counts %>%
      dplyr::mutate(total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
      dplyr::filter(total > min_count)

    # extract txnames
    tx_ids <- counts_filt$TXNAME

    # filter for these transcripts
    bambu_data <- bambu_data[mcols(bambu_data)$transcript_id %in% tx_ids]
  }

  # remove scaffolds
  bambu_data <- bambu_data[grep("chr", seqnames(bambu_data))]

  # remove rows with _ or . in the chromosome
  keep_rows <- !grepl("[._]", seqnames(bambu_data))

  # subset the GRanges
  bambu_data <- bambu_data[keep_rows]

  # filter based on strand
  okstrand <- c("+", "-")
  bambu_data <- bambu_data[strand(bambu_data) %in% okstrand]

  # convert to tibble
  bambu_df <- bambu_data %>% as_tibble()

  # remove version numbers for search
  bambu_df <- bambu_df %>% separate(gene_id, into = "ensg_id", sep = "\\.", remove = FALSE)

  # use mygene to search for gene names
  gene_query <- mygene::queryMany(unique(bambu_df$ensg_id), scopes = "ensembl.gene", fields = "symbol", species = as.character(organism), returnall = TRUE)

  # make df
  gene_df <- as.data.frame(gene_query[["response"]])

  gene_df <- gene_df %>% group_by(query) %>% slice_tail(n = 1) %>% ungroup()

  # if there was no name found, use original ID
  gene_df <- gene_df %>%
    dplyr::mutate(gene_name = case_when(
      is.na(symbol) ~ query,
      !is.na(symbol) ~ symbol
    )) %>%
    dplyr::select(query, gene_name)

  # merge results
  bambu_merged <- merge(bambu_df, gene_df, by.x = "ensg_id", by.y = "query", all.x = T, all.y = F)

  bambu_merged$ensg_id <- NULL

  # make GRanges including new names
  bambu_data_gr <- makeGRangesFromDataFrame(bambu_merged,
                                            keep.extra.columns = TRUE, ignore.strand = FALSE, seqinfo = NULL,
                                            seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand",
                                            starts.in.df.are.0based = FALSE)

  # remove extra mcols
  columns_present <- all(c("gene_name.x", "gene_name.y") %in% colnames(data.frame(bambu_data_gr)))

  if (columns_present) {
    mcols(bambu_data_gr) <- mcols(bambu_data_gr)[, c("source", "type", "score", "phase", "transcript_id", "gene_id", "gene_name.x", "gene_name.y", "exon_number")]
  } else {
    mcols(bambu_data_gr) <- mcols(bambu_data_gr)[, c("source", "type", "score", "phase", "transcript_id", "gene_id", "gene_name", "exon_number")]
  }

  # separate for sorting
  bambu_exons <- bambu_data_gr[bambu_data_gr$type == "exon"]
  bambu_transcripts <- bambu_data_gr[bambu_data_gr$type == "transcript"]

  # sort by chr and locations
  bambu_exons <- sortSeqlevels(bambu_exons)
  bambu_transcripts <- sortSeqlevels(bambu_transcripts)

  bambu_exons <- sort(bambu_exons)
  bambu_transcripts <- sort(bambu_transcripts)

  # recombine for export
  bambu_export <- c(bambu_transcripts, bambu_exons)

  # export filtered gtf
  export(bambu_export, file.path(outdir, "proteome_database_transcripts.gtf"), format = "gtf")

  message("Exported filtered GTF")
}

# function for variant protein sequences
get_variant_orfome <- function(custom_genome, orf_len, txs) {
  custom_genome_fa <- FaFile(custom_genome)
  indexFa(custom_genome)

  # fetch variant transcript sequences
  mut_sequences <- GenomicFeatures::extractTranscriptSeqs(custom_genome_fa, txs)

  # translate ORFs using ORFik
  ORFs <- findMapORFs(txs,
                      mut_sequences,
                      groupByTx = FALSE,
                      longestORF = TRUE,
                      minimumLength = orf_len,
                      startCodon = "ATG",
                      stopCodon = stopDefinition(1))

  ORFs_unlisted <- unlist(ORFs) %>% as_tibble()

  orf_genome_coordinates <- ORFs_unlisted %>%
    rowwise() %>%
    dplyr::mutate(width = end - start) %>%
    group_by(names) %>%
    summarise(chr = seqnames[1],
              start = min(start),
              end = max(end),
              length = sum(width),
              strand = strand[1]) %>%
    ungroup() %>%
    dplyr::select(-length)

  # remove any ORFs from the original ORFs object if they were filtered out due to the length settings above
  ORFs <- ORFs[names(ORFs) %in% orf_genome_coordinates$names]

  orf_seqs <- GenomicFeatures::extractTranscriptSeqs(custom_genome_fa, ORFs)

  # convert nucleotide sequences to amino acid sequences
  orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)

  # create data frame of all possible ORFs
  orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES, protein_sequence = orf_aa_seq, row.names = NULL)

  # remove special characters from the protein sequence and add a column with transcript IDs
  orf_aa_seq_df <- orf_aa_seq_df %>%
    dplyr::mutate(transcript = str_replace(ORF_id, "_.*", "")) %>%
    dplyr::mutate(protein_sequence = str_replace(protein_sequence, "\\*$", ""))

  # combine protein sequences with ORF genomic coordinates
  orf_aa_seq_df_genome_coord <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = c("ORF_id" = "names"))

  fasta_df_mut <- data.frame(
    transcript = names(mut_sequences),
    mut_seq = as.character(mut_sequences),
    stringsAsFactors = FALSE
  )

  return(list(transcript_db = fasta_df_mut, orf_aa_seq_df_genome_coord = orf_aa_seq_df_genome_coord))
}

# fetch variant protein sequences based on VCF file
# function to calculate the number of variable nucleotides
count_variable_nucleotides <- function(mut_seq, wt_seq) {
  sum(strsplit(mut_seq, "")[[1]] != strsplit(wt_seq, "")[[1]])
}

get_variant_protein_seqs <- function(wt_orfome, custom_genome_hm, custom_genome_hm_ht, filteredgtf, genomedb, outdir, orf_len) {
  # import filtered gtf
  txs <- exonsBy(makeTxDbFromGFF(filteredgtf), by = c("tx", "gene"), use.names = TRUE)

  is_homozygous_genome_nonempty <- (file.size(custom_genome_hm) > 0)

  if (is_homozygous_genome_nonempty) {
    orfome_hm <- get_variant_orfome(custom_genome_hm, orf_len, txs)
    orfome_hm_transcript_db <- orfome_hm$transcript_db
    orfome_hm_orf_aa_seq_df_genome_coord <- orfome_hm$orf_aa_seq_df_genome_coord
  }

  orfome_hm_ht <- get_variant_orfome(custom_genome_hm_ht, orf_len, txs)
  orfome_hm_ht_transcript_db <- orfome_hm_ht$transcript_db
  orfome_hm_ht_orf_aa_seq_df_genome_coord <- orfome_hm_ht$orf_aa_seq_df_genome_coord

  # wild-type
  wt_sequences <- GenomicFeatures::extractTranscriptSeqs(genomedb, txs)

  fasta_df_wt <- data.frame(
    transcript = names(wt_sequences),
    wt_seq = as.character(wt_sequences),
    stringsAsFactors = FALSE
  )

  if (is_homozygous_genome_nonempty) {
    variant_transcript_db <- rbind(orfome_hm_transcript_db, orfome_hm_ht_transcript_db)
  } else {
    variant_transcript_db <- orfome_hm_ht_transcript_db
  }

  fasta_df_mut <- distinct(variant_transcript_db)

  # merge wild-type and mutant transcript IDs based on transcript IDs
  fasta_df_merged <- left_join(fasta_df_mut, fasta_df_wt, by = "transcript")

  # apply the count_variable_nucleotides function to the dataframe
  fasta_df_merged$variable_nucleotides <- mapply(count_variable_nucleotides, fasta_df_merged$mut_seq, fasta_df_merged$wt_seq)

  write_tsv(fasta_df_merged, file.path(outdir, "transcriptome_merged.txt"))

  # remove sequences with no variable nucleotides
  fasta_df_merged_mutant <- fasta_df_merged %>%
    filter(variable_nucleotides != 0) %>%
    dplyr::select(transcript, mut_seq) %>%
    unique()

  # subset ORFs for variant transcripts
  orf_aa_seq_df_genome_coord_filtered_ht <- orfome_hm_ht_orf_aa_seq_df_genome_coord %>%
    filter(transcript %in% fasta_df_merged_mutant$transcript) %>%
    dplyr::mutate(mutation_type = "HT")

  if (is_homozygous_genome_nonempty) {
    orf_aa_seq_df_genome_coord_filtered_hm <- orfome_hm_orf_aa_seq_df_genome_coord %>%
      filter(transcript %in% fasta_df_merged_mutant$transcript) %>%
      dplyr::mutate(mutation_type = "HM")

    orf_aa_seq_df_genome_coord_variant <- rbind(orf_aa_seq_df_genome_coord_filtered_hm, orf_aa_seq_df_genome_coord_filtered_ht)
  } else {
    orf_aa_seq_df_genome_coord_variant <- orf_aa_seq_df_genome_coord_filtered_ht
  }

  orf_aa_seq_df_genome_coord_variant <- orf_aa_seq_df_genome_coord_variant %>%
    distinct(protein_sequence, .keep_all = TRUE)

  # add a transcript_id column by splitting ORF_id and taking the first part
  orfome_wt_df <- wt_orfome %>%
    dplyr::mutate(transcript = sapply(strsplit(ORF_id, "_"), `[`, 1)) %>%
    dplyr::select(ORF_sequence, transcript) %>%
    dplyr::mutate(ORF_sequence = str_replace(ORF_sequence, "\\*", ""))

  # remove wild-type proteins
  variant_proteome_flt <- orf_aa_seq_df_genome_coord_variant %>%
    filter(!(protein_sequence %in% orfome_wt_df$ORF_sequence))

  # export protein sequences for the proteome annotation Python script
  variant_protein_seqs <- variant_proteome_flt %>%
    dplyr::mutate(orf_coodinates = paste0(chr, ":", start, "-", end)) %>%
    dplyr::select(transcript, protein_sequence, orf_coodinates, mutation_type) %>%
    unique()

  write_tsv(variant_protein_seqs, file.path(outdir, "Mutant_ORFome_aa.txt"))
}

# generate FASTA of transcript sequences
get_transcript_orfs <- function(filteredgtf, genomedb, orf_len = 30, find_UTR_5_orfs = FALSE, find_UTR_3_orfs = FALSE, referencegtf, outdir) {
  ref_txdb <- NULL
  if ((find_UTR_5_orfs == TRUE) | (find_UTR_5_orfs == TRUE)) {  # required for UTR regions
    ref_txdb <- makeTxDbFromGFF(referencegtf)
  }

  # import filtered gtf
  txs <- exonsBy(makeTxDbFromGFF(filteredgtf), by = c("tx", "gene"), use.names = TRUE)

  # translate into all 3 reading frames
  translate_sequences <- function(sequences) {
    frames <- c(
      sequences,
      subseq(sequences, start = 2),
      subseq(sequences, start = 3)
    )

    sequence_names <- names(sequences)
    aa_names <- c(
      paste0(sequence_names, "_rf1"),
      paste0(sequence_names, "_rf2"),
      paste0(sequence_names, "_rf3")
    )

    aa_sequences <- Biostrings::translate(frames, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    names(aa_sequences) <- aa_names
    return(aa_sequences)
  }

  # ORFik function to find ORFs in GRangesList
  apply_orfik <- function(grl, orfik_min_length, orfik_max_length, orfik_type) {
    # extract transcript nt sequence
    tx_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, grl)
    tx_seq_3nt_indices <- which(width(tx_seqs) >= 3)
    if (length(tx_seq_3nt_indices) == 0) {
      message("None of the transcript sequences contain at least 3 nucleotides, skipping ORF prediction...")
      return(NULL)
    }

    # ORFik
    ORFs <- findMapORFs(grl,
                        tx_seqs,
                        groupByTx = FALSE,
                        longestORF = orfik_type,
                        minimumLength = orfik_min_length,
                        startCodon = "ATG",
                        stopCodon = stopDefinition(1))

    orf_aa_seq_df_genomic_coordinates_merge_frm <- NULL
    ORFs_unlisted <- unlist(ORFs) %>% as_tibble()
    if (NROW(ORFs_unlisted) > 0) {
      # translate all transcript sequences with at least 3 nucleotides into 3 reading frames
      aa_sequences <- translate_sequences(tx_seqs[tx_seq_3nt_indices])

      aa_sequences_df <- data.frame(name = names(aa_sequences), sequence = aa_sequences, stringsAsFactors = FALSE, row.names = NULL) %>%
        mutate(
          tx_id = sub("_.*$", "", name),  # keep up to second dot
          tx_rf_id = sub("^[^.]+\\.[^.]+\\.", "", name)    # keep after second dot
        ) %>% dplyr::select(-name)

      # add width column, filter for ORFs of defined length
      orf_genome_coordinates <- ORFs_unlisted %>%
        rowwise() %>%
        dplyr::mutate(width = end - start) %>%
        group_by(names) %>%
        summarise(chr = seqnames[1],
                  start = min(start),
                  end = max(end),
                  length = sum(width),
                  strand = strand[1]) %>%
        ungroup() %>%
        dplyr::filter(length < (orfik_max_length * 3) - 3) %>% # length ORFs < 30 AA
        dplyr::select(-length)

      # remove any ORFs from original ORF object if they were filtered out due to length settings above
      ORFs <- ORFs[names(ORFs) %in% orf_genome_coordinates$names]

      # rename
      orf_genome_coordinates$ORF_id <- orf_genome_coordinates$names
      orf_genome_coordinates$names <- NULL
      # convert these ORF coordinates into amino acid sequences
      orf_aa_seq <- Biostrings::translate(GenomicFeatures::extractTranscriptSeqs(genomedb, ORFs), if.fuzzy.codon = "solve", no.init.codon = TRUE) # Verify this if it is repetitive
      # create data frame of all possible ORFs
      orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES, ORF_sequence = orf_aa_seq, row.names = NULL)
      # merge with coordinates
      orf_aa_seq_df_genomic_coordinates <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = "ORF_id")
      # separate transcript id
      orf_aa_seq_df_genomic_coordinates <- orf_aa_seq_df_genomic_coordinates %>%
        mutate(ORF_sequence = gsub("\\*", "", ORF_sequence)) %>%
        mutate(tx_id = sub("_\\d+$", "", ORF_id))

      # merge frame dataframe to add reading frame information
      orf_aa_seq_df_genomic_coordinates_merge <- left_join(orf_aa_seq_df_genomic_coordinates, aa_sequences_df, by = "tx_id")

      # get rows containing protein sequences that are sub string of translated sequence to retrieve frame information
      orf_aa_seq_df_genomic_coordinates_merge_frm <- orf_aa_seq_df_genomic_coordinates_merge %>%
        filter(mapply(function(short, long) grepl(short, long, fixed = TRUE), ORF_sequence, sequence))

      orf_aa_seq_df_genomic_coordinates_merge_frm <- orf_aa_seq_df_genomic_coordinates_merge_frm %>%
        mutate(reading_frame = sub(".*_rf", "", tx_rf_id)) %>%
        dplyr::select(-c(sequence, tx_id, tx_rf_id))
    }

    return(orf_aa_seq_df_genomic_coordinates_merge_frm)
  }

  # ORF discovery
  # set ORF max length to large number (to disable)
  # set longestORF to FALSE to ensure we identify CDS
  combined <- apply_orfik(txs, orf_len, 100000, FALSE)

  # extract 5' UTR ORFs
  if (find_UTR_5_orfs == TRUE) {
    # extract 5' UTR regions from ref gtf
    utrs <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    if (length(utrs) > 0) {
      utrs <- utrs[names(utrs) %in% names(txs)]
      if (length(utrs) > 0) {
        # ORF max length is now the main ORF min length
        # only keep longest UTR ORFs
        combined <- rbind(combined, apply_orfik(utrs, 10, orf_len, TRUE))
      } else {
        message("No 5' UTRs from the reference transcript annotation found in the filtered transcripts, skipping ORF prediction in 5' UTRs...")
      }
    } else {
      message("No 5' UTRs found in the reference transcript annotation, skipping ORF prediction in 5' UTRs...")
    }
  }

  # extract 3' UTR ORFs
  if (find_UTR_3_orfs == TRUE) {
    # extract 3' UTR regions from ref gtf
    utrs <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    if (length(utrs) > 0) {
      utrs <- utrs[names(utrs) %in% names(txs)]
      if (length(utrs) > 0) {
        # ORF max length is now the main ORF min length
        # only keep longest UTR ORFs
        combined <- rbind(combined, apply_orfik(utrs, 10, orf_len, TRUE))
      } else {
        message("No 3' UTRs from the reference transcript annotation found in the filtered transcripts, skipping ORF prediction in 3' UTRs...")
      }
    } else {
      message("No 3' UTRs found in the reference transcript annotation, skipping ORF prediction in 3' UTRs...")
    }
  }

  # rename final ORFs with new numerical IDs
  combined <- combined %>%
    separate(ORF_id, into = c("transcript_id"), sep = c("\\_"), remove = F) %>%
    group_by(transcript_id) %>%
    dplyr::mutate(tx_id_number = row_number()) %>%
    ungroup()

  # apply new names
  combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
  combined$transcript_id <- NULL
  combined$tx_id_number <- NULL
  combined <- combined %>% dplyr::mutate(ORF_sequence=str_replace(ORF_sequence, "\\*", ""))

  # export protein seqs for python script
  write_tsv(combined, file.path(outdir, "ORFome_aa.txt"))
  message("Exported ORFik data")

  return(combined)
}

library(optparse)

# define options
option_list = list(
  make_option(c("-g", "--gtf"), type = "character", default = NULL,
              help = "custom user GTF", metavar = "character"),
  make_option(c("-r", "--reference"), type = "character", default = NULL,
              help = "reference GTF", metavar = "character"),
  make_option(c("-G", "--genome"), type = "character", default = NULL,  # optional; only required for generating variant-aware databases
              help = "reference genome FASTA", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "transcript counts file", metavar = "character"),
  make_option(c("-m", "--mincount"), type = "integer", default = NULL,
              help = "Minimum transcript count", metavar = "integer"),
  make_option(c("-o", "--organism"), type = "character", default = NULL,
              help = "Organism", metavar = "character"),
  make_option(c("-l", "--length"), type = "integer", default = NULL,
              help = "Minimum ORF length", metavar = "integer"),
  make_option(c("-u", "--uorfs"), type = "logical", default = NULL,
              help = "Find uORFs", metavar = "TRUE/FALSE"),
  make_option(c("-d", "--dorfs"), type = "logical", default = NULL,
              help = "Find dORFs", metavar = "TRUE/FALSE"),
  make_option(c("-v", "--vcf"), type = "character", default = NULL,
              help = "VCF file", metavar = "character"),
  make_option(c("-s", "--savepath"), type = "character", default = NULL,
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# define inputs
gtf_path <- opt$gtf
reference_gtf <- opt$reference
tx_count_path <- opt$counts
minimum_tx_count <- opt$m
ref_genome <- opt$genome
organism <- opt$organism
min_orf_length <- opt$length
find_5_orfs <- opt$uorfs
find_3_orfs <- opt$dorfs
vcf_file <- opt$vcf
output_directory <- opt$savepath

# check if any required arguments are missing

error_message <- ""

if (is.null(gtf_path)) {
  error_message <- "Please provide a custom user GTF."
} else if (is.null(reference_gtf)) {
  error_message <- "Please provide a reference GTF."
} else if (is.null(organism)) {
  error_message <- "Please specify one of the following organisms (case-insensitive): 'HUMAN', 'MOUSE', 'RAT', 'CAEEL' (C. elegans), 'DROME' (Drosophila), or 'DANRE' (Zebrafish)."
} else if (is.null(min_orf_length)) {
  error_message <- "Please specify the minimum ORF length."
} else if (is.null(find_5_orfs)) {
  error_message <- "Please indicate whether ORFs from 5' UTRs should be predicted."
} else if (is.null(find_3_orfs)) {
  error_message <- "Please indicate whether ORFs from 3' UTRs should be predicted."
} else if (is.null(output_directory)) {
  error_message <- "Please specify an output directory."
}

if (nchar(error_message) > 0) {
  stop(error_message)
}

# check the input arguments

organism <- toupper(trimws(as.character(organism)))
min_orf_length <- floor(as.numeric(min_orf_length))
if (!is.null(minimum_tx_count) && !is.null(tx_count_path)) {
  minimum_tx_count <- floor(as.numeric(minimum_tx_count))
}

if (!file.exists(gtf_path)) {
  error_message <- paste0("The custom user GTF '", gtf_path, "' does not exist.")
} else if (!file.exists(reference_gtf)) {
  error_message <- paste0("The reference GTF '", reference_gtf, "' does not exist.")
} else if (!is.null(ref_genome) && !file.exists(ref_genome)) {
  error_message <- paste0("The reference genome FASTA '", ref_genome, "' does not exist.")
} else if (!is.null(tx_count_path) && !file.exists(tx_count_path)) {
  error_message <- paste0("The transcript counts file '", tx_count_path, "' does not exist.")
} else if (!is.null(vcf_file) && !file.exists(vcf_file)) {
  error_message <- paste0("The VCF file '", vcf_file, "' does not exist.")
} else if (!(organism %in% c("HUMAN", "MOUSE", "CAEEL", "DROME", "RAT", "DANRE"))) {
  error_message <- paste0("Organism '", as.character(opt$organism), "' is not supported. Please specify one of the following (case-insensitive): 'HUMAN', 'MOUSE', 'RAT', 'CAEEL' (C. elegans), 'DROME' (Drosophila), or 'DANRE' (Zebrafish).")
} else if (file.exists(output_directory) && !dir.exists(output_directory)) {
  error_message <- paste0("'", output_directory, "' exists but is not a directory.")
} else if (!(find_5_orfs %in% c(TRUE, FALSE))) {
  error_message <- "Please specify 'TRUE' or 'FALSE' on whether to predict ORFs from 5' UTRs."
} else if (!(find_3_orfs %in% c(TRUE, FALSE))) {
  error_message <- "Please specify 'TRUE' or 'FALSE' on whether to predict ORFs from 3' UTRs."
} else if (!(is.finite(min_orf_length) && min_orf_length >= 0)) {
  error_message <- paste0("The minimum ORF length '", opt$length, "' cannot be parsed as a non-negative integer.")
} else if (!is.null(minimum_tx_count) && !is.null(tx_count_path) && !(is.finite(minimum_tx_count) && minimum_tx_count >= 0)) {
  error_message <- paste0("The minimum transcript count '", opt$m, "' cannot be parsed as a non-negative integer.")
}

if (nchar(error_message) > 0) {
  stop(error_message)
}

# select the reference genome corresponding to the specified organism
if (organism == "HUMAN") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  genomedb <- BSgenome.Hsapiens.UCSC.hg38
} else if (organism == "MOUSE") {
  library(BSgenome.Mmusculus.UCSC.mm39)
  genomedb <- BSgenome.Mmusculus.UCSC.mm39
} else if (organism == "CAEEL") {   # C. elegans
  library(BSgenome.Celegans.UCSC.ce11)
  genomedb <- BSgenome.Celegans.UCSC.ce11
} else if (organism == "DROME") {   # Drosophila
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
} else if (organism == "RAT") {
  library(BSgenome.Rnorvegicus.UCSC.rn7)
  genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
} else if (organism == "DANRE") {   # Zebrafish
  library(BSgenome.Drerio.UCSC.danRer11)
  genomedb <- BSgenome.Drerio.UCSC.danRer11
}

# create the output directory if it does not already exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# load the rest of the libraries needed for the remainder of the script to run
source("global.R")

if (!grepl("proteome_database_transcripts.gtf", gtf_path)) {
  if (!is.null(tx_count_path)) {
    filter_custom_gtf(customgtf = gtf_path, organism = organism, tx_counts = tx_count_path, min_count = minimum_tx_count, outdir = output_directory)
  } else {
    filter_custom_gtf(customgtf = gtf_path, organism = organism, outdir = output_directory)
  }
  filtered_gtf <- file.path(output_directory, "proteome_database_transcripts.gtf")
} else {
  filtered_gtf <- gtf_path
}

wt_orfome <- get_transcript_orfs(filteredgtf = filtered_gtf, genomedb = genomedb, orf_len = min_orf_length, find_UTR_5_orfs = find_5_orfs, find_UTR_3_orfs = find_3_orfs, referencegtf = reference_gtf, outdir = output_directory)

# generate custom genome
if (!is.null(ref_genome) && !is.null(vcf_file)) { # if a genome and VCF have been uploaded
  # define bash script command to inject variants into genome
  custom_genome_command <- paste0("bash bin/database_module/generate_custom_genome.sh -g ", ref_genome, " -r ", reference_gtf, " -v ", vcf_file, " -o ", output_directory)
  system(custom_genome_command)

  # fetch variant protein sequences based on variants provided in VCF file
  genome_alt_hm    <- file.path(output_directory, basename(ref_genome) %>% str_replace(., ".fa",    "_hm.fa"))
  genome_alt_hm_ht <- file.path(output_directory, basename(ref_genome) %>% str_replace(., ".fa", "_hm_ht.fa"))

  if (!file.exists(genome_alt_hm)) {
    message(paste0("The homozygous consensus sequence '", genome_alt_hm, "' does not exist. Skip finding variant protein sequences..."))
  } else if (!file.exists(genome_alt_hm_ht)) {
    message(paste0("The heterozygous consensus sequence '", genome_alt_hm_ht, "' does not exist. Skip finding variant protein sequences..."))
  } else if (file.size(genome_alt_hm_ht) == 0) {
    message(paste0("The heterozygous consensus sequence '", genome_alt_hm_ht, "' is empty. Skip finding variant protein sequences..."))
  } else {
    # apply variant protein function
    get_variant_protein_seqs(wt_orfome = wt_orfome, custom_genome_hm = genome_alt_hm, custom_genome_hm_ht = genome_alt_hm_ht, filteredgtf = file.path(output_directory, "proteome_database_transcripts.gtf"), genomedb = genomedb, outdir = output_directory, min_orf_length)
  }
}