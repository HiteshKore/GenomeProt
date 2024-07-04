# functions

# import gtf and filter for minimum transcript counts
filter_custom_gtf <- function(customgtf, tx_counts=NA, min_count=NA) {
  
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
    counts_filt <- counts %>% dplyr::filter(if_any(where(is.numeric), ~ . >= min_count))
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
  
  # export filtered gtf
  export(bambu_data, "ORFome_transcripts.gtf", format="gtf")

}

# export FASTA of transcript sequences
get_transcript_seqs <- function (filteredgtf, organism) {
  
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
  
  # export
  writeXStringSet(tx_seqs, "ORFome_transcripts_nt.fasta")
  
}

# import proteomics file (fragpipe, maxquant etc)
import_proteomics_data <- function(proteomics_file) {
  
  # required columns
  # Peptide, Protein, Mapped Proteins
  # optional: Protein Start
  
  #proteomics_file <- "~/Documents/proteogenomics/2024/ben_data/report.pr_matrix.tsv"
  protfile <- fread(proteomics_file)
  
  if (c("Protein Start") %in% colnames(protfile)) {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Protein Start`, `Mapped Proteins`)
  } else {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Mapped Proteins`)
  }
  
  # add a column of every mapped ORF
  protfile <- protfile %>% mutate(all_mappings = case_when(
    (is.na(`Mapped Proteins`) == TRUE | nchar(`Mapped Proteins`) == 0) ~ paste0(Protein),
    TRUE ~ paste0(Protein, ", ", `Mapped Proteins`))
  )
  
  # separate into one row per mapped ORF
  prot_expanded <- separate_rows(protfile, all_mappings, sep = "\\, |\\;")
  prot_expanded <- prot_expanded %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = fixed(" "))) %>% 
    dplyr::mutate(PID = str_replace(all_mappings, "\\,", "\\."))
  prot_expanded <- prot_expanded[!(base::duplicated(prot_expanded)),]
  
  if (c("Protein Start") %in% colnames(prot_expanded)) {
    prot_expanded <- prot_expanded %>% dplyr::select(PID, Peptide, `Protein Start`)
  } else {
    prot_expanded <- prot_expanded %>% dplyr::select(PID, Peptide)
  }
  
  # remove reverse mappings?
  # should we report on reverse mappings?
  
  prot_expanded <- prot_expanded %>% 
    dplyr::filter(!startsWith(PID, "rev"))
  
  prot_expanded$peptide <- prot_expanded$Peptide
  prot_expanded$Peptide <- NULL
  
  return(prot_expanded)
  
}

# import metadata (old)
import_metadata <- function(metadata_file, proteomics_data) {
  
  # required columns
  # accession (protein), orf_genomic_coordinates, orf_transcript_coordinates, transcript?
  
  #metadata_file <- "~/Documents/proteogenomics/2024/miguel_tx_and_proteomics/ProteomeDb_metadata.txt"
  #proteomics_data <- pd
  
  # import database of ORF ids
  metadata <- fread(metadata_file)
  metadata$accession <- gsub("\\,", "\\.", metadata$accession)
  metadata$PID <- paste0(metadata$accession, "|CO=", metadata$orf_genomic_coordinates)
  metadata <- separate(metadata, orf_genomic_coordinates, into = c("chromosome", "start", "end"), sep = "[:-]", remove = FALSE)
  metadata <- separate(metadata, orf_transcript_coordinates, into = c("txstart", "txend"), sep = "[-]", remove = TRUE)
  
  metadata <- metadata %>% 
    mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
  
  metadata <- metadata %>% 
    dplyr::filter(txstart != 0)
  
  # merge database with proteomics data
  df <- merge(proteomics_data, metadata, by="PID", all.x=TRUE, all.y=FALSE)
  
  df <- df %>% dplyr::filter(!is.na(transcript))
  
  df$orf_tx_id <- paste0(df$PID, "|", df$transcript)
  
  # make lower case
  df$peptide <- df$Peptide
  df$Peptide <- NULL
  
  # find peptide sequence exact matches in protein sequences
  find_exact_match <- function(peptide, protein_sequence) {
    match_positions <- gregexpr(pattern = peptide, text = protein_sequence)
    start_position <- match_positions[[1]][1]
    if (start_position == -1) {
      return(-1)
    } else {
      return(start_position)
    }
  }
  
  # find leftover peptides with no exact match, sequences now with one mismatch in protein sequences
  find_one_mismatch <- function(peptide, protein_sequence) {
    peptide_length <- nchar(peptide)
    protein_length <- nchar(protein_sequence)
    
    for (i in 1:(protein_length - peptide_length + 1)) {
      sub_seq <- substr(protein_sequence, i, i + peptide_length - 1)
      dist <- stringdist::stringdist(peptide, sub_seq, method = "hamming")
      if (dist == 1) {
        return(i)
      }
    }
    return(-1)
  }
  
  # combine both functions
  find_peptide_position <- function(peptide, protein_sequence) {
    exact_match_position <- find_exact_match(peptide, protein_sequence)
    
    if (exact_match_position != -1) {
      return(exact_match_position)
    } else {
      return(find_one_mismatch(peptide, protein_sequence))
    }
  }
  
  # apply to dataframe
  df$calc_pep_start <- mapply(find_peptide_position, df$peptide, df$protein_sequence)
  
  # if exact match or one mismatch isn't found, use proteomics defined start site
  if (c("Protein Start") %in% colnames(df)) {
    df <- df %>% mutate(pep_start = case_when(
      (calc_pep_start != `Protein Start` & calc_pep_start != -1) ~ calc_pep_start,
      (calc_pep_start != `Protein Start` & calc_pep_start == -1) ~ `Protein Start`,
      TRUE ~ calc_pep_start
    ))
  } else {
    df$pep_start <- df$calc_pep_start
  }
  
  # get peptide end location within every ORF
  df$pep_end <- df$pep_start + nchar(df$peptide)
  df$calc_pep_start <- NULL
  
  return(df)
  
}

# import custom database FASTA of ORFs
import_fasta <- function(fasta_file, proteomics_data, gtf_file) {
  
  # required columns
  # accession (protein), orf_genomic_coordinates, orf_transcript_coordinates, transcript?
  #fasta_file <- fasta_import_file
  #fasta_file <- "~/Documents/proteogenomics/2024/miguel_tx_and_proteomics/fragpipe/ProteomeDb.fasta"
  #proteomics_data <- pd
  
  db <- readAAStringSet(fasta_file, format="fasta", use.names=TRUE)
  
  # convert to df
  convert_AA_to_df <- function(AAstring) {
    data.frame(name=names(AAstring),
               protein_sequence=as.character(AAstring),
               row.names=NULL)
  }
  
  # convert
  df <- convert_AA_to_df(db)
  
  # NOTE: FASTA headers belonging to multiple genomic locations are removed...
  df <- separate(df, name, into = c("PID"), sep = "\\ |\\,", remove = FALSE)
  df <- separate(df, name, into = c("header", "gene_id", "transcript_id", "strand"), sep = "\\ ", remove = FALSE)
  df <- separate(df, gene_id, into = c("gene_id"), sep = "\\,", remove = FALSE)
  df <- separate(df, strand, into = c("strand"), sep = "\\,", remove = FALSE)
  
  df <- separate(df, PID, into = c("protein_name", "location"), sep = "\\|CO=", remove = FALSE)
  df <- separate(df, location, into = c("chromosome", "start", "end"), sep = "\\:|-", remove = FALSE)
  
  df <- df %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("GN="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("TA="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("ST="))) 
  
  df_expanded <- separate_rows(df, transcript_id, sep = "\\,|\\, ")
  df_expanded <- df_expanded[!(base::duplicated(df_expanded)),]
  
  df_filtered <- df_expanded %>% dplyr::filter(PID %in% proteomics_data$PID)
  
  # Get df of orf genomic coords
  orf_genomic_coords_df <- df_filtered %>% dplyr::select(PID, protein_name, chromosome, start, end, strand, transcript_id)
  
  # create unique id
  orf_genomic_coords_df$uniq_id <- paste0(orf_genomic_coords_df$protein_name, "_", orf_genomic_coords_df$transcript_id)
  
  # split into list based on ORF and transcript id
  # can't just use mapToTranscripts, as some transcript have multiple orfs
  orf_genomic_coords_list <- split(orf_genomic_coords_df, orf_genomic_coords_df$uniq_id)
  
  # get transcripts for mapping
  txs <- transcripts(gtf_file)
  names(txs) <- txs$tx_name
  
  # NOTE
  # this is super slow
  # maybe just do separately for transcripts with multiple orfs?
  # do pmapToTranscripts for transcripts with 1 orf
  
  create_granges_and_map_to_txs <- function(df_list_obj) {
    #df_list_obj <- orf_genomic_coords_list[[1]]
    genomic_coords <- df_list_obj %>% dplyr::filter(!is.na(chromosome)) %>% makeGRangesFromDataFrame(
      keep.extra.columns=TRUE,
      ignore.strand=FALSE,
      seqinfo=NULL,
      seqnames.field="chromosome",
      start.field="start",
      end.field="end",
      strand.field="strand",
      starts.in.df.are.0based=FALSE,
      na.rm=FALSE)
    
    txs_filtered <- txs[names(txs) %in% genomic_coords$transcript_id] # filter for only transcripts with peptides
    transcript_coords <- mapToTranscripts(genomic_coords, txs_filtered)
    transcript_coords_df <- transcript_coords %>% as_tibble()
    
    transcript_coords_df <- transcript_coords_df %>% dplyr::select(seqnames, start, end)
    
    transcript_coords_df$transcript_id <- transcript_coords_df$seqnames
    transcript_coords_df$seqnames <- NULL
    
    transcript_coords_df$PID <- df_list_obj$PID
    transcript_coords_df$txstart <- transcript_coords_df$start
    transcript_coords_df$txend <- transcript_coords_df$end
    transcript_coords_df$end <- NULL
    transcript_coords_df$start <- NULL
    
    return(transcript_coords_df)
  }
  
  orf_transcriptomic_coords_list <- lapply(orf_genomic_coords_list[1:10], create_granges_and_map_to_txs)
  
  orf_transcriptomic_coords <- do.call("rbind", orf_transcriptomic_coords_list)
  
  ## --- ##
  
  # orf_transcriptomic_coords <- orf_transcriptomic_coords %>% 
  #   dplyr::filter(start >= 1)
  
  # merge database with proteomics data
  df_merged <- merge(df_filtered, orf_transcriptomic_coords, by=c("PID", "transcript_id"), all.x=T, all.y=F)
  df_merged <- df_merged[!(base::duplicated(df_merged)),]
  
  df_merged <- df_merged %>% dplyr::filter(!is.na(transcript_id))
  
  metadata <- merge(df_merged, proteomics_data, by="PID", all.x=F, all.y=T)
  metadata <- metadata[!(base::duplicated(metadata)),]
  
  metadata <- metadata %>% dplyr::filter(!is.na(transcript_id))
  
  # find peptide sequence exact matches in protein sequences
  find_exact_match <- function(peptide, protein_sequence) {
    match_positions <- gregexpr(pattern = peptide, text = protein_sequence)
    start_position <- match_positions[[1]][1]
    if (start_position == -1) {
      return(-1)
    } else {
      return(start_position)
    }
  }
  
  # find leftover peptides with no exact match, sequences now with one mismatch in protein sequences
  find_one_mismatch <- function(peptide, protein_sequence) {
    peptide_length <- nchar(peptide)
    protein_length <- nchar(protein_sequence)
    
    for (i in 1:(protein_length - peptide_length + 1)) {
      sub_seq <- substr(protein_sequence, i, i + peptide_length - 1)
      dist <- stringdist::stringdist(peptide, sub_seq, method = "hamming")
      if (dist == 1) {
        return(i)
      }
    }
    return(-1)
  }
  
  # combine both functions
  find_peptide_position <- function(peptide, protein_sequence) {
    exact_match_position <- find_exact_match(peptide, protein_sequence)
    
    if (exact_match_position != -1) {
      return(exact_match_position)
    } else {
      return(find_one_mismatch(peptide, protein_sequence))
    }
  }
  
  # apply to dataframe
  metadata$calc_pep_start <- mapply(find_peptide_position, metadata$peptide, metadata$protein_sequence)
  
  # if exact match or one mismatch isn't found, use proteomics defined start site
  if (c("Protein Start") %in% colnames(metadata)) {
    metadata <- metadata %>% mutate(pep_start = case_when(
      (calc_pep_start != `Protein Start` & calc_pep_start != -1) ~ calc_pep_start,
      (calc_pep_start != `Protein Start` & calc_pep_start == -1) ~ `Protein Start`,
      TRUE ~ calc_pep_start
    ))
  } else {
    metadata$pep_start <- metadata$calc_pep_start
  }
  
  # get peptide end location within every ORF
  metadata$pep_end <- metadata$pep_start + nchar(metadata$peptide)
  metadata$calc_pep_start <- NULL
  
  return(metadata)
  
}

# get peptide coords from transcript coords
extract_peptide_coords <- function(metadata_df, tx_orfs) {
  
  # subset metadata for conversion of peptide coords
  subset_df <- metadata_df %>% dplyr::select(orf_tx_id, transcript, PID, gene, pep_start, pep_end, peptide)
  txcoordsdf_subset <- merge(subset_df, tx_orfs, by=c("orf_tx_id", "transcript", "gene"), all.x=F, all.y=F)
  txcoordsdf_subset <- txcoordsdf_subset[!(base::duplicated(txcoordsdf_subset)),]
  
  # get peptide locations within every transcript
  txcoordsdf_subset$txstart <- as.numeric(txcoordsdf_subset$txstart)
  txcoordsdf_subset$txend <- as.numeric(txcoordsdf_subset$txend)
  
  txcoordsdf_subset$start <- txcoordsdf_subset$txstart + ((txcoordsdf_subset$pep_start-1) * 3)
  txcoordsdf_subset$end <- txcoordsdf_subset$start + (nchar(txcoordsdf_subset$peptide) * 3)
  txcoordsdf_subset$width <- txcoordsdf_subset$end - txcoordsdf_subset$start
  txcoordsdf_subset$peptide <- as.factor(txcoordsdf_subset$peptide)
  
  # GRanges from df of transcriptomic coords of peptides
  coords_peptides <- makeGRangesFromDataFrame(txcoordsdf_subset,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field="transcript",
                                              start.field="start",
                                              end.field="end",
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE,
                                              na.rm=TRUE)
  # set names
  names(coords_peptides) <- c(coords_peptides$transcript)
  
  return(coords_peptides)
  
}






