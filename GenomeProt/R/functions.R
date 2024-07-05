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
    #counts_filt <- counts %>% dplyr::filter(if_any(where(is.numeric), ~ . >= min_count))
    counts_filt <- counts %>% 
      mutate(total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
      dplyr::filter(total > min_count)
    
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

# import custom database FASTA of ORFs
import_fasta <- function(fasta_file, proteomics_data, gtf_file) {
  
  # fasta_file <- fasta_import_file
  # gtf_file <- gtf
  # proteomics_data <- pd
  
  # get transcripts for mapping
  txs <- exonsBy(gtf_file, by=c("tx"), use.names=TRUE)
  
  # import fasta
  db <- readAAStringSet(fasta_file, format="fasta", use.names=TRUE)
  
  # convert to df
  convert_AA_to_df <- function(AAstring) {
    data.frame(name=names(AAstring),
               protein_sequence=as.character(AAstring),
               row.names=NULL)
  }
  
  # convert
  df <- convert_AA_to_df(db)
  
  # NOTE: FASTA headers belonging to multiple genomic locations are removed
  # format headers into columns
  df <- df %>%
    separate(name, into = c("PID"), sep = "\\ |\\,", remove = FALSE) %>%
    separate(name, into = c("header", "gene_id", "transcript_id", "strand"), sep = "\\ ", remove = FALSE) %>%
    separate(gene_id, into = c("gene_id"), sep = "\\,", remove = FALSE) %>%
    separate(strand, into = c("strand"), sep = "\\,", remove = FALSE) %>%
    separate(PID, into = c("protein_name", "location"), sep = "\\|CO=", remove = FALSE) %>%
    separate(location, into = c("chromosome", "start", "end"), sep = "\\:|-", remove = FALSE)
  
  df <- df %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("GN="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("TA="))) %>% 
    dplyr::mutate(across(where(is.character), str_remove_all, pattern = c("ST="))) 
  
  # add one row per transcript
  df_expanded <- separate_rows(df, transcript_id, sep = "\\,|\\, ")
  df_expanded <- df_expanded[!(base::duplicated(df_expanded)),]
  
  # filter for proteins with mapped peptides
  df_expanded <- df_expanded %>% dplyr::filter(PID %in% proteomics_data$PID)
  
  # get df of orf genomic coords
  orf_genomic_coords_df <- df_expanded %>% dplyr::select(PID, protein_name, chromosome, start, end, strand, transcript_id)
  
  # list per transcript
  orf_genomic_coords_list <- split(orf_genomic_coords_df, orf_genomic_coords_df$transcript_id)
  
  create_granges_and_map_to_txs <- function(df_list_obj) {
    #df_list_obj <- orf_genomic_coords_list["ENSMUST00000025707.9"][[1]]
    genomic_coords <- df_list_obj %>% dplyr::filter(!is.na(chromosome)) %>% makeGRangesFromDataFrame(
      keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
      seqnames.field="chromosome", start.field="start", end.field="end",
      strand.field="strand", starts.in.df.are.0based=FALSE, na.rm=FALSE)
    
    txs_filtered <- txs[names(txs) %in% genomic_coords$transcript_id] # filter for singe transcript
    
    transcript_coords <- pmapToTranscripts(genomic_coords, txs_filtered) # map to transcript coords using exons
    
    mcols(transcript_coords)$PID <- mcols(genomic_coords)$PID # add back PID
    mcols(transcript_coords)$protein_name <- mcols(genomic_coords)$protein_name # add back protein_name
    
    transcript_coords_df <- transcript_coords %>%
      as_tibble() %>%
      mutate(
        transcript_id = seqnames,
        txstart = start,
        txend = end,
        orf_length = width
      ) %>%
      dplyr::select(-seqnames, -start, -end, -strand, -width)
      
    return(transcript_coords_df)
      
  }
  
  orf_transcriptomic_coords_list <- lapply(orf_genomic_coords_list, create_granges_and_map_to_txs)
  
  orf_transcriptomic_coords <- do.call("rbind", orf_transcriptomic_coords_list)
  
  # ORFs that could not be mapped to transcripts are returned as 0 length
  orf_transcriptomic_coords <- orf_transcriptomic_coords %>% 
    dplyr::filter(txstart > 1 & orf_length > 0)
  
  ## --- ##
  
  # merge database with proteomics data
  df_merged <- merge(df_expanded, orf_transcriptomic_coords, by=c("PID", "transcript_id", "protein_name"), all.x=T, all.y=F)
  df_merged <- df_merged[!(base::duplicated(df_merged)),]
  
  df_merged <- df_merged %>% dplyr::filter(!is.na(transcript_id))
  
  metadata <- merge(df_merged, proteomics_data, by="PID", all.x=F, all.y=T)
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
  metadata$mapped_pep_start <- NULL
  
  return(metadata)
  
}

# get peptide coords from transcript coords
extract_peptide_coords <- function(metadata_df, tx_orfs) {

  # subset metadata for conversion of peptide coords
  subset_df <- metadata_df %>% dplyr::select(orf_tx_id, transcript_id, PID, gene_id, pep_start, pep_end, peptide)
  txcoordsdf_subset <- merge(subset_df, tx_orfs, by=c("orf_tx_id", "transcript_id", "gene_id"), all.x=F, all.y=F)
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
                                              seqnames.field="transcript_id",
                                              start.field="start",
                                              end.field="end",
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE,
                                              na.rm=TRUE)
  # set names
  names(coords_peptides) <- c(coords_peptides$transcript)
  
  return(coords_peptides)
  
}






