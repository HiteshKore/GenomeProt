
# define functions

source("global.R")

# import gtf and filter for minimum transcript counts
filter_custom_gtf <- function(customgtf, tx_counts=NA, min_count=NA) {
  # customgtf <- "~/Documents/GenomeProt_tmp/test_datasets/db_module/miguel_subset.gtf"
  # tx_counts <- "~/Documents/GenomeProt_tmp/test_datasets/db_module/miguel_counts.csv"
  # min_count <- 40
  
  # import bambu gtf
  bambu_data <- rtracklayer::import(customgtf)
  
  if (!missing(tx_counts)) {
    
    # define default value
    if (missing(min_count)) {
      min_count <- 5
    }
    
    # read in counts
    counts <- fread(tx_counts)
    
    # filter for counts greater than or equal to min_count in any numeric column
    counts_filt <- counts %>% 
      mutate(total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
      dplyr::filter(total > as.numeric(min_count))
    
    # extract txnames
    tx_ids <- counts_filt$TXNAME
    # filter for these transcripts
    bambu_data <- bambu_data[mcols(bambu_data)$transcript_id %in% tx_ids]
    
  } 
  
  # remove scaffolds etc
  bambu_data <- bambu_data[grep("chr", seqnames(bambu_data))]
  
  # filter based on strand
  okstrand <- c("+", "-")
  bambu_data <- bambu_data[strand(bambu_data) %in% okstrand]
  
  # remove extra mcols
  mcols(bambu_data) <- mcols(bambu_data)[, c("source", "type", "score", "phase", "transcript_id", "gene_id", "exon_number")]
  bambu_exons <- bambu_data[bambu_data$type == "exon"]
  bambu_transcripts <- bambu_data[bambu_data$type == "transcript"]
  
  # sort by chr and locations
  bambu_exons <- sortSeqlevels(bambu_exons)
  bambu_transcripts <- sortSeqlevels(bambu_transcripts)
  
  bambu_exons <- sort(bambu_exons)
  bambu_transcripts <- sort(bambu_transcripts)
  
  bambu_export <- c(bambu_transcripts, bambu_exons)
  
  # export filtered gtf
  export(bambu_export, "database_output/proteome_database_transcripts.gtf", format="gtf")
  
  print("Exported filtered GTF")
  
}

# export FASTA of transcript sequences
get_transcript_orfs <- function (filteredgtf, organism, orf_len=30, find_UTR_5_orfs=FALSE, find_UTR_3_orfs=FALSE, referencegtf) {
  # filteredgtf <- "~/Documents/GenomeProt_tmp/demo_datasets/3_db_module/bambu_transcript_annotations.gtf"
  # organism <- "human"
  # orf_len <- 30
  # find_UTR_5_orfs <- FALSE
  # find_UTR_3_orfs <- TRUE
  # referencegtf <- "~/Documents/gencode_annotations/gencode.v44.annotation.gtf"
  
  # import filtered gtf as a txdb
  txdb <- makeTxDbFromGFF(filteredgtf)
  txs <- exonsBy(txdb, by=c("tx","gene"), use.names=TRUE)
  
  # convert to GRangesList
  txs.granges <- GRangesList(txs)
  
  # use this GTF to extract transcript sequences based on genome (extracts novel sequences too)
  if (organism=="mouse") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39::Mmusculus
    tx_seqs <- extractTranscriptSeqs(genomedb, txs.granges)
  } else if (organism=="human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
    tx_seqs <- extractTranscriptSeqs(genomedb, txs.granges)
  }
  
  # get nt sequence as df
  #fasta_export_nucleotides <- data.frame(transcript = names(tx_seqs), sequence = tx_seqs, row.names=NULL)
  #writeXStringSet(tx_seqs, "database_output/ORFome_transcripts_nt.fasta")
  
  # ORFik to find ORFs
  ORFs <- findMapORFs(txs.granges, 
                      tx_seqs, 
                      groupByTx = FALSE, 
                      longestORF = TRUE, 
                      minimumLength = as.numeric(orf_len), 
                      startCodon = "ATG",
                      stopCodon = stopDefinition(1))
  
  # unlist GRL
  ORFs_unlisted <- unlist(ORFs) %>% as_tibble()
  
  orf_genome_coordinates <- ORFs_unlisted %>% group_by(names) %>% 
    summarise(chr = seqnames[1],
              start = min(start),
              end = max(end),
              strand = strand[1])
  
  orf_genome_coordinates$ORF_id <- orf_genome_coordinates$names
  
  # convert these ORFs coordinates into nucleotide sequences
  orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, ORFs)
  
  # convert the nucleotide sequences to amino acid sequences
  orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
  
  # create data frame of all possible ORFs
  orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,ORF_sequence = orf_aa_seq, row.names=NULL) 
  
  orf_aa_seq_df_genomic_coordinates <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = "ORF_id")
  orf_aa_seq_df_genomic_coordinates$names <- NULL
  
  combined <- orf_aa_seq_df_genomic_coordinates

  if (find_UTR_5_orfs == FALSE & find_UTR_3_orfs == FALSE) {
    
    ref_transcripts <- rtracklayer::import(referencegtf)
    export(ref_transcripts, "database_output/ref_transcripts_in_data.gtf", format="gtf")
    
  } else if (find_UTR_5_orfs == TRUE & find_UTR_3_orfs == FALSE) {
    
    ref_transcripts <- rtracklayer::import(referencegtf)
    ref_transcripts <- ref_transcripts[mcols(ref_transcripts)$transcript_id %in% names(txs)]
    
    # need to add if statement, what if no reference transcripts are found in the data?
    export(ref_transcripts, "database_output/ref_transcripts_in_data.gtf", format="gtf")
    
    ref_txdb <- makeTxDbFromGFF("database_output/ref_transcripts_in_data.gtf")
    
    utrs5 <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs5_filtered <- utrs5[names(utrs5) %in% names(txs)]
    
    # extract UTR transcript sequences
    utr_seqs <- extractTranscriptSeqs(genomedb, utrs5_filtered)
    
    # ORFik to find ORFs
    utrORFs <- findMapORFs(utrs5_filtered,
                           utr_seqs, 
                           groupByTx = FALSE, 
                           longestORF = TRUE,
                           minimumLength = 10,
                           startCodon = "ATG",
                           stopCodon = stopDefinition(1))
    
    utrORFs_unlisted <- unlist(utrORFs) %>% as_tibble()
    
    utrORFs_unlisted$new_names <- paste0(utrORFs_unlisted$names, "_utr")
    
    utr_orf_genome_coordinates <- utrORFs_unlisted %>% 
      rowwise() %>% 
      dplyr::mutate(width = end - start) %>% 
      group_by(new_names, names) %>% 
      summarise(chr = seqnames[1],
                start = min(start),
                end = max(end),
                length = sum(width),
                strand = strand[1]) %>% 
      ungroup() %>% 
      dplyr::filter(length < (as.numeric(orf_len)*3)-3) %>% # length ORFs < 30 AA
      dplyr::select(-length)
    
    utrORFs <- utrORFs[names(utrORFs) %in% utr_orf_genome_coordinates$names]
    
    utr_orf_genome_coordinates$ORF_id <- utr_orf_genome_coordinates$names
    
    # convert these ORFs coordinates into nucleotide sequences
    utr_orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, utrORFs)
    
    # convert the nucleotide sequences to amino acid sequences
    utr_orf_aa_seq <- Biostrings::translate(utr_orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    utr_orf_aa_seq_df <- data.frame(ORF_id = utr_orf_aa_seq@ranges@NAMES, 
                                    ORF_sequence = utr_orf_aa_seq, row.names=NULL)
    
    utr_orf_aa_seq_df_genomic_coordinates <- left_join(utr_orf_aa_seq_df, utr_orf_genome_coordinates, by = "ORF_id")
    
    utr_orf_aa_seq_df_genomic_coordinates$ORF_id <- utr_orf_aa_seq_df_genomic_coordinates$new_names
    utr_orf_aa_seq_df_genomic_coordinates$new_names <- NULL
    utr_orf_aa_seq_df_genomic_coordinates$names <- NULL
    
    combined <- NULL
    
    combined <- rbind(orf_aa_seq_df_genomic_coordinates, utr_orf_aa_seq_df_genomic_coordinates)
    
    combined <- combined %>% 
      separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
      group_by(transcript_id) %>% 
      mutate(tx_id_number = row_number()) %>% 
      ungroup()
    
    combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
    
    combined$transcript_id <- NULL
    combined$tx_id_number <- NULL
    
  } else if (find_UTR_5_orfs == FALSE & find_UTR_3_orfs == TRUE) {
    
    ref_transcripts <- rtracklayer::import(referencegtf)
    ref_transcripts <- ref_transcripts[mcols(ref_transcripts)$transcript_id %in% names(txs)]
    
    # need to add if statement, what if no reference transcripts are found in the data?
    export(ref_transcripts, "database_output/ref_transcripts_in_data.gtf", format="gtf")
    
    ref_txdb <- makeTxDbFromGFF("database_output/ref_transcripts_in_data.gtf")
    
    utrs3 <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs3_filtered <- utrs3[names(utrs3) %in% names(txs)]
    
    # extract UTR transcript sequences
    utr_seqs <- extractTranscriptSeqs(genomedb, utrs3_filtered)
    
    # ORFik to find ORFs
    utrORFs <- findMapORFs(utrs3_filtered,
                           utr_seqs,
                           groupByTx = FALSE, 
                           longestORF = TRUE,
                           minimumLength = 10,
                           startCodon = "ATG",
                           stopCodon = stopDefinition(1))
    
    utrORFs_unlisted <- unlist(utrORFs) %>% as_tibble()
    
    utrORFs_unlisted$new_names <- paste0(utrORFs_unlisted$names, "_utr")
    
    utr_orf_genome_coordinates <- utrORFs_unlisted %>% 
      rowwise() %>% 
      dplyr::mutate(width = end - start) %>% 
      group_by(new_names, names) %>% 
      summarise(chr = seqnames[1],
                start = min(start),
                end = max(end),
                length = sum(width),
                strand = strand[1]) %>% 
      ungroup() %>% 
      dplyr::filter(length < (as.numeric(orf_len)*3)-3) %>% # length ORFs < 30 AA
      dplyr::select(-length)
    
    utrORFs <- utrORFs[names(utrORFs) %in% utr_orf_genome_coordinates$names]
    
    utr_orf_genome_coordinates$ORF_id <- utr_orf_genome_coordinates$names
    
    # convert these ORFs coordinates into nucleotide sequences
    utr_orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, utrORFs)
    
    # convert the nucleotide sequences to amino acid sequences
    utr_orf_aa_seq <- Biostrings::translate(utr_orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    utr_orf_aa_seq_df <- data.frame(ORF_id = utr_orf_aa_seq@ranges@NAMES, 
                                    ORF_sequence = utr_orf_aa_seq, row.names=NULL)
    
    utr_orf_aa_seq_df_genomic_coordinates <- left_join(utr_orf_aa_seq_df, utr_orf_genome_coordinates, by = "ORF_id")
    
    utr_orf_aa_seq_df_genomic_coordinates$ORF_id <- utr_orf_aa_seq_df_genomic_coordinates$new_names
    utr_orf_aa_seq_df_genomic_coordinates$new_names <- NULL
    utr_orf_aa_seq_df_genomic_coordinates$names <- NULL
    
    combined <- NULL
    
    combined <- rbind(orf_aa_seq_df_genomic_coordinates, utr_orf_aa_seq_df_genomic_coordinates)
    
    combined <- combined %>% 
      separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
      group_by(transcript_id) %>% 
      mutate(tx_id_number = row_number()) %>% 
      ungroup()
    
    combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
    
    combined$transcript_id <- NULL
    combined$tx_id_number <- NULL
    
  } else if (find_UTR_5_orfs == TRUE & find_UTR_3_orfs == TRUE) {
    
    ref_transcripts <- rtracklayer::import(referencegtf)
    ref_transcripts <- ref_transcripts[mcols(ref_transcripts)$transcript_id %in% names(txs)]
    
    # need to add if statement, what if no reference transcripts are found in the data?
    export(ref_transcripts, "database_output/ref_transcripts_in_data.gtf", format="gtf")
    
    ref_txdb <- makeTxDbFromGFF("database_output/ref_transcripts_in_data.gtf")
    
    utrs5 <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs3 <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs_both <- c(utrs5, utrs3)
    
    utrsboth_filtered <- utrs_both[names(utrs_both) %in% names(txs)]
    
    # extract UTR transcript sequences
    utr_seqs <- extractTranscriptSeqs(genomedb, utrsboth_filtered)
    
    # ORFik to find ORFs
    utrORFs <- findMapORFs(utrsboth_filtered,
                           utr_seqs,
                           groupByTx = FALSE, 
                           longestORF = TRUE,
                           minimumLength = 10,
                           startCodon = "ATG",
                           stopCodon = stopDefinition(1))
    
    utrORFs_unlisted <- unlist(utrORFs) %>% as_tibble()
    
    utrORFs_unlisted$new_names <- paste0(utrORFs_unlisted$names, "_utr")
    
    utr_orf_genome_coordinates <- utrORFs_unlisted %>% 
      rowwise() %>% 
      dplyr::mutate(width = end - start) %>% 
      group_by(new_names, names) %>% 
      summarise(chr = seqnames[1],
                start = min(start),
                end = max(end),
                length = sum(width),
                strand = strand[1]) %>% 
      ungroup() %>% 
      dplyr::filter(length < (as.numeric(orf_len)*3)-3) %>% # length ORFs < 30 AA
      dplyr::select(-length)
    
    utrORFs <- utrORFs[names(utrORFs) %in% utr_orf_genome_coordinates$names]
    
    utr_orf_genome_coordinates$ORF_id <- utr_orf_genome_coordinates$names
    
    # convert these ORFs coordinates into nucleotide sequences
    utr_orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, utrORFs)
    
    # convert the nucleotide sequences to amino acid sequences
    utr_orf_aa_seq <- Biostrings::translate(utr_orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    utr_orf_aa_seq_df <- data.frame(ORF_id = utr_orf_aa_seq@ranges@NAMES, 
                                    ORF_sequence = utr_orf_aa_seq, row.names=NULL)
    
    utr_orf_aa_seq_df_genomic_coordinates <- left_join(utr_orf_aa_seq_df, utr_orf_genome_coordinates, by = "ORF_id")
    
    utr_orf_aa_seq_df_genomic_coordinates$ORF_id <- utr_orf_aa_seq_df_genomic_coordinates$new_names
    utr_orf_aa_seq_df_genomic_coordinates$new_names <- NULL
    utr_orf_aa_seq_df_genomic_coordinates$names <- NULL
    
    combined <- NULL
    
    combined <- rbind(orf_aa_seq_df_genomic_coordinates, utr_orf_aa_seq_df_genomic_coordinates)
    
    combined <- combined %>% 
      separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
      group_by(transcript_id) %>% 
      mutate(tx_id_number = row_number()) %>% 
      ungroup()
    
    combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
    
    combined$transcript_id <- NULL
    combined$tx_id_number <- NULL
  }
  
  # export protein seqs
  write_tsv(combined, "database_output/ORFome_aa.txt")
  
  print("Exported ORFik data")
  
}

option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="custom user GTF", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="reference user GTF", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL,
              help="transcript counts file", metavar="character"),
  make_option(c("-m", "--mincount"), type="integer", default=NULL,
              help="Minimum transcript count", metavar="integer"),
  make_option(c("-o", "--organism"), type="character", default=NULL,
              help="Organism", metavar="character"),
  make_option(c("-l", "--length"), type="integer", default=NULL,
              help="Minimum ORF length", metavar="integer"),
  make_option(c("-u", "--uorfs"), type="logical", default=NULL,
              help="Find uORFs", metavar="TRUE/FALSE"),
  make_option(c("-d", "--dorfs"), type="logical", default=NULL,
              help="Find dORFs", metavar="TRUE/FALSE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

gtf_path <- opt$gtf
reference_gtf <- opt$reference
tx_count_path <- opt$counts
minimum_tx_count <- opt$m
organism <- opt$organism
min_orf_length <- opt$length
find_5_orfs <- opt$uorfs
find_3_orfs <- opt$dorfs

# make output directory
system("mkdir database_output")

# run filter_custom_gtf, check if counts are present
if (!is.null(tx_count_path)) {
  filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, tx_counts=tx_count_path, min_count=minimum_tx_count)
} else {
  filtered_gtf <- filter_custom_gtf(customgtf=gtf_path)
}

# run orfik to find ORFs
get_transcript_orfs(filteredgtf="database_output/proteome_database_transcripts.gtf", organism=organism, orf_len=min_orf_length, find_UTR_5_orfs=find_5_orfs, find_UTR_3_orfs=find_3_orfs, referencegtf=reference_gtf)



