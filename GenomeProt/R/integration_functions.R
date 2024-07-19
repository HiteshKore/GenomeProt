# functions

# import proteomics file (fragpipe, maxquant etc)
import_proteomics_data <- function(proteomics_file) {
  
  # required columns
  # Peptide, Protein, Mapped Proteins
  # optional: Protein Start
  
  # if its metamorpheus:
  # Peptide column is called 'Base Sequence', proteins matched are 'Protein Groups'
  # # Use str_replace to capture groups and replace with a delimiter
  # # delimits based on TA= (any num chars) followed by a |
  # split_data <- str_replace_all(data, "(TA=[A-Za-z0-9.]+)\\|", "\\1_SPLIT_")
  # 
  # # Split based on the new delimiter
  # split_data <- str_split(split_data, "_SPLIT_", simplify = TRUE)
  # 
  # # Convert to a data frame
  # df <- as.data.frame(split_data, stringsAsFactors = FALSE)
  
  
  #proteomics_file <- "~/Documents/proteogenomics/2024/ben_data/report.pr_matrix.tsv"
  protfile <- fread(proteomics_file)
  
  if (c("Protein Start") %in% colnames(protfile)) {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Protein Start`, `Mapped Proteins`)
  } else {
    protfile <- protfile %>% dplyr::select(Peptide, Protein, `Mapped Proteins`)
  }
  
  # add a column of every mapped ORF
  protfile <- protfile %>% dplyr::mutate(
    all_mappings = case_when(
      (is.na(`Mapped Proteins`) == TRUE | nchar(`Mapped Proteins`) == 0) ~ paste0(Protein),
      TRUE ~ paste0(Protein, ", ", `Mapped Proteins`)
      )
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

# import custom database FASTA of ORFs
import_fasta <- function(fasta_file, proteomics_data, gtf_file) {
  
  # fasta_file <- fasta_import_file
  # gtf_file <- gtf
  # proteomics_data <- pd
  
  tx_lengths <- transcriptLengths(gtf_file)
  tx_lengths$transcript_id <- tx_lengths$tx_name
  tx_lengths <- tx_lengths %>% dplyr::select(transcript_id, tx_len)
  
  # get transcripts for mapping
  txs <- exonsBy(gtf_file, by=c("tx"), use.names=T)
  
  # import fasta
  db <- readAAStringSet(fasta_file, format="fasta", use.names=TRUE)
  
  # convert to df
  convert_AA_to_df <- function(AAstring) {
    data.frame(name=names(AAstring),
               protein_sequence=as.character(AAstring),
               protein_length=as.numeric(nchar(as.character(AAstring))),
               row.names=NULL)
  }
  
  # convert
  df <- convert_AA_to_df(db)
  df$nt_length <- df$protein_length * 3
  
  # format headers into columns
  suppressWarnings({
    df <- df %>%
      separate(name, into = c("PID"), sep = "\\ |\\,", remove = FALSE) 
    
    # filter for proteins with mapped peptides
    df <- df %>% dplyr::filter(PID %in% proteomics_data$PID)
    
    df <- df %>% 
      separate(name, into = c("header", "gene_id", "transcript_id", "strand"), sep = "\\ ", remove = FALSE) %>%
      separate(gene_id, into = c("gene_id"), sep = "\\,", remove = FALSE) %>%
      separate(strand, into = c("strand"), sep = "\\,", remove = FALSE) %>%
      separate(PID, into = c("protein_name", "location"), sep = "\\|CO=", remove = FALSE) %>%
      separate(location, into = c("chromosome", "start", "end"), sep = "\\:|-", remove = FALSE)
  })
  
  df <- df %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("GN="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("TA="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("ST="))) 
  
  # add one row per transcript
  df_expanded <- separate_rows(df, transcript_id, sep = "\\,|\\, ")
  df_expanded <- df_expanded[!(base::duplicated(df_expanded)),]
  
  df_expanded <- merge(df_expanded, tx_lengths, by="transcript_id", all.x=T, all.y=F)
  
  df_expanded$start <- as.numeric(df_expanded$start)
  df_expanded$end <- as.numeric(df_expanded$end)
  
  # get df of orf genomic coords
  orf_genomic_coords_df <- df_expanded
  
  orf_genomic_coords_df <- orf_genomic_coords_df %>% 
    rowwise() %>% 
    mutate(
      stranded_start = case_when(
        strand == "+" ~ min(start,end),
        strand == "-" ~ max(start,end),
      ))
  
  txs_filtered <- txs[names(txs) %in% orf_genomic_coords_df$transcript_id]
  txs_unlisted <- unlist(txs_filtered)
  
  # function to convert genomic start position to transcript position
  convert_gene_pos_to_transcript_pos <- function(txdb, transcript_id, genomic_start, strand) {
    
    # Get the exons for the specific transcript
    # transcript_id <- "ENST00000585111.2"
    # genomic_start <- 64506119
    # strand <- "-"

    tx_exons <- txdb[names(txdb) == transcript_id]
    
    if (length(tx_exons) == 0) {
      transcript_start <- -1
    }
    
    # Find the cumulative transcript position for each exon
    exon_data <- data.frame(
      exon_start = start(tx_exons),
      exon_end = end(tx_exons),
      exon_rank = mcols(tx_exons)$exon_rank
    ) %>%
      arrange(exon_rank) %>%
      mutate(
        exon_length = exon_end - exon_start + 1,
        total_length = cumsum(exon_length)
      )
    
    # Find which exon the genomic start position falls into
    exon_index <- which(genomic_start >= exon_data$exon_start & genomic_start <= exon_data$exon_end)
    
    if (length(exon_index) == 0) {
      transcript_start <- -1
    } else {
    
      # calculate the transcript start position
      if (strand == "+") {
        if (exon_index == 1) {
          transcript_start <- genomic_start - exon_data$exon_start[exon_index]
        } else {
          transcript_start <- (genomic_start - exon_data$exon_start[exon_index]) + exon_data$total_length[exon_index - 1]
        }
      } else {
        if (exon_index == 1) {
          transcript_start <- exon_data$exon_end[exon_index] - genomic_start
        } else {
          # edit...
          transcript_start <- (exon_data$exon_end[exon_index] - genomic_start) + exon_data$total_length[exon_index - 1]
        }
      }
    }
    return(transcript_start)
  }
  
  # wrapper function to apply to df
  apply_conversion_to_df <- function(txdb, df) {
    df <- df %>%
      rowwise() %>%
      mutate(txstart = convert_gene_pos_to_transcript_pos(txdb, transcript_id, stranded_start, strand))
    
    return(df)
  }
  
  # apply to df
  #orf_genomic_coords_df_subset <- orf_genomic_coords_df[1:5000,]
  orf_transcriptomic_coords <- apply_conversion_to_df(txs_unlisted, orf_genomic_coords_df)
  
  # ORFs that could not be mapped to transcripts are returned with -1 starts
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txstart >= 1)
  # remove ORF coords outside transcript ends
  orf_transcriptomic_coords$txend <- orf_transcriptomic_coords$txstart + orf_transcriptomic_coords$nt_length - 1
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% dplyr::filter(txend < tx_len)
  
  # combine with proteomics data
  
  metadata <- merge(orf_transcriptomic_coords, proteomics_data, by="PID", all.x=T, all.y=T)
  metadata <- metadata[!(base::duplicated(metadata)),]
  
  metadata <- metadata %>% dplyr::filter(!is.na(transcript_id) & !is.na(txstart))
  
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
  metadata$mapped_pep_start <- mapply(find_peptide_position, metadata$peptide, metadata$protein_sequence)
  
  # if exact match or one mismatch isn't found, use proteomics defined start site
  if (c("Protein Start") %in% colnames(metadata)) {
    metadata <- metadata %>% mutate(pep_start = case_when(
      (mapped_pep_start != `Protein Start` & mapped_pep_start != -1) ~ mapped_pep_start,
      (mapped_pep_start != `Protein Start` & mapped_pep_start == -1) ~ `Protein Start`,
      TRUE ~ mapped_pep_start
    ))
  } else {
    metadata$pep_start <- metadata$mapped_pep_start
  }
  
  # get peptide end location within every ORF
  metadata$pep_end <- metadata$pep_start + nchar(metadata$peptide)
  
  metadata <- metadata %>% dplyr::filter(pep_start < protein_length & pep_end < protein_length)
  
  metadata$mapped_pep_start <- NULL
  
  return(metadata)
  
}

# get peptide coords from transcript coords
extract_peptide_coords <- function(metadata_df) {

  # subset metadata for conversion of peptide coords
  txcoordsdf_subset <- metadata_df %>% dplyr::select(orf_tx_id, transcript_id, PID, gene_id, pep_start, pep_end, peptide, txstart, txend)
  
  # get peptide locations within every transcript
  txcoordsdf_subset$txstart <- as.numeric(txcoordsdf_subset$txstart)
  txcoordsdf_subset$txend <- as.numeric(txcoordsdf_subset$txend)
  
  txcoordsdf_subset$start <- txcoordsdf_subset$txstart + ((txcoordsdf_subset$pep_start-1) * 3) - 1
  txcoordsdf_subset$end <- txcoordsdf_subset$start + (nchar(txcoordsdf_subset$peptide) * 3)
  txcoordsdf_subset$width <- txcoordsdf_subset$end - txcoordsdf_subset$start
  txcoordsdf_subset$peptide <- as.factor(txcoordsdf_subset$peptide)
  
  txcoordsdf_subset <- txcoordsdf_subset %>% 
    dplyr::filter(start > txstart & end < txend)
  
  # GRanges from df of transcriptomic coords of peptides
  coords_peptides <- makeGRangesFromDataFrame(txcoordsdf_subset,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field="transcript_id",
                                              start.field="start",
                                              end.field="end",
                                              strand.field="strand",
                                              starts.in.df.are.0based=TRUE,
                                              na.rm=TRUE)
  # set names
  names(coords_peptides) <- c(coords_peptides$transcript)
  
  
  # filter peptides for bad mappings
  coords_peptides <- subset(coords_peptides, start(coords_peptides) != 0)
  coords_peptides <- subset(coords_peptides, end(coords_peptides) < coords_peptides$txend)
  coords_peptides <- subset(coords_peptides, (mcols(coords_peptides)$txstart != mcols(coords_peptides)$pep_start))
  coords_peptides <- subset(coords_peptides, (mcols(coords_peptides)$txend != mcols(coords_peptides)$pep_end))
  
  return(coords_peptides)
  
}






