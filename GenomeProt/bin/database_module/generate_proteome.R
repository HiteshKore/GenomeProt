source("global.R")

# import gtf and filter for minimum transcript counts
filter_custom_gtf <- function(customgtf, organism, tx_counts=NA, min_count=NA, outdir) {
  
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
      dplyr::filter(total > as.numeric(min_count))
    
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
  bambu_df <- bambu_df %>% separate(gene_id, into="ensg_id", sep="\\.", remove = FALSE)
  
  # use mygene to search for gene names
  gene_query <- mygene::queryMany(unique(bambu_df$ensg_id), scopes="ensembl.gene", fields="symbol", species=as.character(organism),  returnall=TRUE)
  
  # make df
  gene_df <- as.data.frame(gene_query[["response"]])
  
  gene_df<-gene_df %>% group_by(query) %>% slice_tail(n=1)%>%ungroup()
  
  # if there was no name found, use original ID
  gene_df <- gene_df %>% 
    dplyr::mutate(gene_name = case_when(
      is.na(symbol) ~ query,
      !is.na(symbol) ~ symbol
    )) %>% 
    dplyr::select(query, gene_name)
  
  # merge results
  bambu_merged <- merge(bambu_df, gene_df, by.x="ensg_id", by.y="query", all.x=T, all.y=F)
  
  bambu_merged$ensg_id <- NULL
  
  # make GRanges including new names
  bambu_data_gr <- makeGRangesFromDataFrame(bambu_merged,
                                            keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                            seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand",
                                            starts.in.df.are.0based=FALSE)
  
  # remove extra mcols
  columns_present <- all(c("gene_name.x", "gene_name.y") %in%  colnames(data.frame(bambu_data_gr)))
  
  if (columns_present) {
    mcols(bambu_data_gr) <- mcols(bambu_data_gr)[, c("source", "type", "score", "phase", "transcript_id", "gene_id", "gene_name.x","gene_name.y", "exon_number")]
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
  export(bambu_export, paste0(outdir, "/proteome_database_transcripts.gtf"), format="gtf")
  
  print("Exported filtered GTF")
  
}

#function for variant protein seqs
get_variant_orfome<-function(custom_genome,custom_gtf,orf_len,txs_grl){
  custom_genome_hm_fa <- FaFile(custom_genome)
  indexFa(custom_genome)
  
  # fetch variant transcript sequences
  mut_sequences <- GenomicFeatures::extractTranscriptSeqs(custom_genome_hm_fa, txs_grl)
  
  # translate ORFs using ORFik
  ORFs <- findMapORFs(txs_grl,
                      mut_sequences, 
                      groupByTx = FALSE,
                      longestORF = TRUE, 
                      minimumLength = as.numeric(orf_len), 
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
  # remove any ORFs from original ORF object if they were filtered out due to length settings above
  ORFs <- ORFs[names(ORFs) %in% orf_genome_coordinates$names]
  
  
  
  
  orf_seqs <- GenomicFeatures::extractTranscriptSeqs(custom_genome_hm_fa, ORFs)
  
  
  # convert nucleotide sequences to amino acid sequences
  orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
  
  # create data frame of all possible ORFs
  orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,protein_sequence = orf_aa_seq, row.names=NULL)
  
  
  
  # remove special characters from protein sequence and add column with transcript id
  orf_aa_seq_df <- orf_aa_seq_df %>%
    dplyr::mutate(transcript=str_replace(ORF_id, "_.*", "")) %>%
    dplyr::mutate(protein_sequence=str_replace(protein_sequence,"\\*$", ""))
  
  # combine protein sequences with ORF genomic coordinates
  orf_aa_seq_df_genome_coord <- left_join(orf_aa_seq_df,orf_genome_coordinates,by=c("ORF_id"="names"))
  # 
  fasta_df_mut <- data.frame(
    transcript = names(mut_sequences),
    mut_seq = as.character(mut_sequences),
    stringsAsFactors = FALSE
  )
  return(list(transcript_db=fasta_df_mut,orf_aa_seq_df_genome_coord=orf_aa_seq_df_genome_coord))
}

# fetch variant protein sequences based on VCF file
get_variant_protein_seqs <- function(wt_orfome, custom_genome_hm, custom_genome_hm_ht,custom_gtf, organism,outdir, orf_len) {
  
  
  if (file.exists(custom_genome_hm) && file.exists(custom_genome_hm_ht)) {
    
    # import filtered gtf as a txdb
    txdb <- makeTxDbFromGFF(custom_gtf)
    txs <- exonsBy(txdb, by=c("tx", "gene"), use.names=TRUE)
    # convert txdb to GRangesList
    txs_grl <- GRangesList(txs)
    
    orfome_hm<-get_variant_orfome(custom_genome_hm,custom_gtf,orf_len,txs_grl)
    orfome_hm_transcript_db<-orfome_hm$transcript_db
    orfome_hm_orf_aa_seq_df_genome_coord<-orfome_hm$orf_aa_seq_df_genome_coord
    
    orfome_hm_ht<-get_variant_orfome(custom_genome_hm_ht,custom_gtf,orf_len,txs_grl)
    orfome_hm_ht_transcript_db<-orfome_hm_ht$transcript_db
    orfome_hm_ht_orf_aa_seq_df_genome_coord<-orfome_hm_ht$orf_aa_seq_df_genome_coord
    
    
    variant_proteome <- rbind(orfome_hm_orf_aa_seq_df_genome_coord,orfome_hm_ht_orf_aa_seq_df_genome_coord)
    variant_proteome <- distinct(variant_proteome)
    
    #rbind variant transcripts
    variant_transcript_db<-rbind(orfome_hm_transcript_db,orfome_hm_ht_transcript_db)
    fasta_df_mut <- distinct(variant_transcript_db)
    
    #write output for annotation script
    #write_tsv(variant_proteome, paste0(outdir, "variant_proteome.txt"))
    
    # wild type
    # set organism
    if (organism == "HUMAN") {
      library(BSgenome.Hsapiens.UCSC.hg38)
      genomedb <- BSgenome.Hsapiens.UCSC.hg38
    } else if (organism == "MOUSE") {
      library(BSgenome.Mmusculus.UCSC.mm39)
      genomedb <- BSgenome.Mmusculus.UCSC.mm39
    } else if (organism == "CAEEL") { #celegans
      library(BSgenome.Celegans.UCSC.ce11)
      genomedb <- BSgenome.Celegans.UCSC.ce11
    } else if (organism == "DROME") { #drosophila
      library(BSgenome.Dmelanogaster.UCSC.dm6)
      genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
    } else if (organism == "RAT") { #rat
      library(BSgenome.Rnorvegicus.UCSC.rn7)
      genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
    } else if (organism == "DANRE") { #zebrafish
      library(BSgenome.Drerio.UCSC.danRer11)
      genomedb <- BSgenome.Drerio.UCSC.danRer11
    }
    
    wt_genomedb <- genomedb
    
    
    wt_sequences <- extractTranscriptSeqs(wt_genomedb, txs_grl)
    
    fasta_df_wt <- data.frame(
      transcript = names(wt_sequences),
      wt_seq = as.character(wt_sequences),
      stringsAsFactors = FALSE
    )
    
    
    # merge wild type and mutant transcript ids based on transcript ids
    fasta_df_merged <- left_join(fasta_df_mut,fasta_df_wt,by="transcript")
    
    # # function to calculate the number of variable nucleiotides
    count_variable_nucleotides <- function(mut_seq, wt_seq) {
      sum(strsplit(mut_seq, "")[[1]] != strsplit(wt_seq, "")[[1]])
    }
    # 
    # # apply the function to the dataframe
    fasta_df_merged$variable_nucleotides <- mapply(count_variable_nucleotides, fasta_df_merged$mut_seq, fasta_df_merged$wt_seq)
    
    write_tsv(fasta_df_merged, paste0(outdir, "transcriptome_merged.txt"))
    
    # # remove sequences with no variable nucleotides
    fasta_df_merged_mutant <- fasta_df_merged %>%
      filter(variable_nucleotides != 0) %>% 
      dplyr::select(transcript,mut_seq) %>%
      unique()
    
    # subset ORFs for variant transcripts
    orf_aa_seq_df_genome_coord_filtered_hm <- orfome_hm_orf_aa_seq_df_genome_coord %>%
      filter(transcript %in% fasta_df_merged_mutant$transcript) %>%dplyr::mutate(mutation_type="HM")
    
    orf_aa_seq_df_genome_coord_filtered_ht <- orfome_hm_ht_orf_aa_seq_df_genome_coord %>%
      filter(transcript %in% fasta_df_merged_mutant$transcript) %>%dplyr::mutate(mutation_type="HT")
    
    orf_aa_seq_df_genome_coord_variant<-rbind(orf_aa_seq_df_genome_coord_filtered_hm,orf_aa_seq_df_genome_coord_filtered_ht)%>%distinct(protein_sequence,.keep_all = TRUE)
    
    
    
    
    
    # adding the transcript_id column by splitting ORF_id and taking the first part
    orfome_wt_df <- wt_orfome %>% 
      dplyr::mutate(transcript = sapply(strsplit(ORF_id, "_"), `[`, 1)) %>%
      dplyr::select(ORF_sequence,transcript)%>%dplyr::mutate(ORF_sequence=str_replace(ORF_sequence,"\\*", ""))
    #remove wild type proteins
    variant_proteome_flt <- orf_aa_seq_df_genome_coord_variant %>%
      filter(!(protein_sequence %in% orfome_wt_df$ORF_sequence))
    
    #format variant_proteome_flt
    variant_protein_seqs <- variant_proteome_flt %>%
      dplyr::mutate(orf_coodinates=paste0(chr,":",start,"-",end)) %>%
      dplyr::select(transcript,protein_sequence,orf_coodinates,mutation_type) %>%
      unique()
    
    # export protein seqs for python script
    write_tsv(variant_protein_seqs, paste0(outdir, "/Mutant_ORFome_aa.txt"))
    
  } else {
    
    print("File does not exist")
    
  }
  
}

# generate FASTA of transcript sequences
get_transcript_orfs <- function(filteredgtf, organism, orf_len=30, find_UTR_5_orfs=FALSE, find_UTR_3_orfs=FALSE, referencegtf, outdir) {
  
  # set organism
  if (organism == "HUMAN") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
  } else if (organism == "MOUSE") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39
  } else if (organism == "CAEEL") {
    library(BSgenome.Celegans.UCSC.ce11)
    genomedb <- BSgenome.Celegans.UCSC.ce11
  } else if (organism == "DROME") {
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
  } else if (organism == "RAT") {
    library(BSgenome.Rnorvegicus.UCSC.rn7)
    genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
  } else if (organism == "DANRE") {
    library(BSgenome.Drerio.UCSC.danRer11)
    genomedb <- BSgenome.Drerio.UCSC.danRer11
  }
  
  # required for UTR regions
  ref_txdb <- makeTxDbFromGFF(referencegtf)
  
  # import filtered gtf as a txdb
  txdb <- makeTxDbFromGFF(filteredgtf)
  txs <- exonsBy(txdb, by=c("tx", "gene"), use.names=TRUE)
  
  # convert txdb to GRangesList
  txs_grl <- GRangesList(txs)
  
  
  ######
  # translate into all 3 reading frames 
  translate_sequence <- function(sequence, transcript_name) {
    sequence<-Biostrings::DNAString(sequence)
    if( nchar(sequence)<3) return(NULL)
    
    frames <- list(
      subseq(sequence, start = 1),
      subseq(sequence, start = 2),
      subseq(sequence, start = 3)
    )
    
    aa_sequences <- lapply(seq_along(frames), function(i) {
      
      aa <-as.character(Biostrings::translate(frames[[i]], if.fuzzy.codon = "solve", no.init.codon = TRUE))
      names(aa) <- paste0(transcript_name, "_rf", as.character(i))
      return(aa)
    })
    
    return(aa_sequences)
  }
  
  
  # ORFik function to find ORFs in GRangesList
  apply_orfik <- function(grl, orfik_min_length, orfik_max_length, orfik_type) {
    orf_aa_seq_df_genomic_coordinates<-NULL
    # extract transcript nt sequence
    tx_seqs <- extractTranscriptSeqs(genomedb, grl)
    
    # # ORFik
    ORFs <- findMapORFs(grl,
                        tx_seqs, 
                        groupByTx = FALSE,
                        longestORF = orfik_type, 
                        minimumLength = as.numeric(orfik_min_length), 
                        startCodon = "ATG",
                        stopCodon = stopDefinition(1))
    ORFs_unlisted <- unlist(ORFs) %>% as_tibble()
    
    
    if(NROW(ORFs_unlisted) >0){
      
      # translate the transcripts sequences into 3 reading frame to extract the frame information
      
      aa_sequences_list <- lapply(seq_along(tx_seqs), function(i) {
        translate_sequence(tx_seqs[[i]], names(tx_seqs)[i])
      })
      
      
      aa_sequences <- unlist(aa_sequences_list)
      
      aa_sequences_df <- data.frame(name = names(aa_sequences), sequence = aa_sequences, stringsAsFactors = FALSE, row.names = NULL) %>%
        mutate(
          tx_id = sub("_.*$", "", name),  # keep up to second dot
          tx_rf_id = sub("^[^.]+\\.[^.]+\\.", "", name)    # keep after second dot
        )%>%dplyr::select(-name)
      
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
        dplyr::filter(length < (as.numeric(orfik_max_length)*3)-3) %>% # length ORFs < 30 AA
        dplyr::select(-length)
      
      # remove any ORFs from original ORF object if they were filtered out due to length settings above
      ORFs <- ORFs[names(ORFs) %in% orf_genome_coordinates$names]
      
      # rename
      orf_genome_coordinates$ORF_id <- orf_genome_coordinates$names
      orf_genome_coordinates$names <- NULL
      # convert these ORF coordinates into nucleotide sequences
      orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, ORFs)
      #convert the nucleotide sequences to amino acid sequences
      orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
      # create data frame of all possible ORFs
      orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,ORF_sequence = orf_aa_seq, row.names=NULL) 
      # merge with coordinates
      orf_aa_seq_df_genomic_coordinates <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = "ORF_id")
      #separate transcript id
      orf_aa_seq_df_genomic_coordinates<-orf_aa_seq_df_genomic_coordinates%>%mutate(ORF_sequence=gsub("\\*","",ORF_sequence))%>%
        mutate(tx_id=sub("_\\d+$","",ORF_id))
      #merge frame dataframe to add reading frame information
      orf_aa_seq_df_genomic_coordinates_merge<-left_join(orf_aa_seq_df_genomic_coordinates,aa_sequences_df,by="tx_id")
      #get rows containing protein sequences that are sub string of translated sequence to retrieve frame information
      orf_aa_seq_df_genomic_coordinates_merge_frm<-orf_aa_seq_df_genomic_coordinates_merge%>%filter( mapply(function(short,long)grepl(short,long,fixed=TRUE),ORF_sequence,sequence) )
      
      # match_vec <- vapply(seq_len(nrow(orf_aa_seq_df_genomic_coordinates_merge)), function(i) {
      #   grepl(orf_aa_seq_df_genomic_coordinates_merge$ORF_sequence[i],
      #         orf_aa_seq_df_genomic_coordinates_merge$sequence[i],
      #         fixed = TRUE)
      # }, logical(1))
      # 
      # orf_aa_seq_df_genomic_coordinates_merge_frm <- orf_aa_seq_df_genomic_coordinates_merge[match_vec, ] %>%
      #   mutate(reading_frame = sub(".*_rf", "", tx_rf_id)) %>%
      #   dplyr::select(-c(sequence, tx_id, tx_rf_id))
      #print(orf_aa_seq_df_genomic_coordinates_merge_frm)
      print(dim(orf_aa_seq_df_genomic_coordinates_merge_frm))
      orf_aa_seq_df_genomic_coordinates_merge_frm<-orf_aa_seq_df_genomic_coordinates_merge_frm%>%mutate(reading_frame=sub(".*_rf","",tx_rf_id))%>%dplyr::select(-c(sequence,tx_id,tx_rf_id))
      
      
      
      
    }
    
    
    return(orf_aa_seq_df_genomic_coordinates_merge_frm)
    
  }
  # ORF discovery
  # set ORF max length to large number (to disable)
  # set longestORF to FALSE to ensure we identify CDS
  all_ORFs <- apply_orfik(txs_grl, as.numeric(orf_len), 100000, FALSE)
  #print(all_ORFs)
  
  
  # create tmp copy
  tmp <- all_ORFs
  # extract 5UTR ORFs
  if (find_UTR_5_orfs == TRUE) {
    
    # extract 5' UTR regions from ref gtf
    utrs5 <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    
    if (length(utrs5) > 0){
      utrs5_filtered <- utrs5[names(utrs5) %in% names(txs)]
      
      # ORF max length is now the main ORF min length
      # only keep longest UTR ORFs
      utr5_ORFs <- apply_orfik(utrs5_filtered, 10, as.numeric(orf_len), TRUE)
      # update tmp
      tmp <- rbind(tmp, utr5_ORFs)
      
    }
    
    
  }
  # extract 3UTR ORFs
  if (find_UTR_3_orfs == TRUE) {
    
    # extract 3' UTR regions from ref gtf
    utrs3 <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    if (length(utrs3) > 0){
      utrs3_filtered <- utrs3[names(utrs3) %in% names(txs)]
      
      # ORF max length is now the main ORF min length
      # only keep longest UTR ORFs
      utr3_ORFs <- apply_orfik(utrs3_filtered, 10, as.numeric(orf_len), TRUE)
      
      # update tmp
      tmp <- rbind(tmp, utr3_ORFs)
    }
    
  }
  
  
  # 
  # # rename final ORFs with new numerical IDs
  combined <- tmp %>% 
    separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
    group_by(transcript_id) %>% 
    dplyr::mutate(tx_id_number = row_number()) %>% 
    ungroup()
  
  # # apply new names
  combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
  combined$transcript_id <- NULL
  combined$tx_id_number <- NULL
  
  combined<-combined%>%dplyr::mutate(ORF_sequence=str_replace(ORF_sequence,"\\*", ""))
  
  
  
  # export protein seqs for python script
  write_tsv(combined, paste0(outdir, "/ORFome_aa.txt"))
  print("Exported ORFik data")
  
  return(combined)
  
}

# define options
option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="custom user GTF", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="reference GTF", metavar="character"),
  make_option(c("-G", "--genome"), type="character", default=NULL, #Optional.Only requires for generating variant aware database
              help="reference genome fa", metavar="character"),
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
              help="Find dORFs", metavar="TRUE/FALSE"),
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="vcf file", metavar="character"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
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



# Check if the path ends with a "/"
if (substr(output_directory, nchar(output_directory), nchar(output_directory)) != "/") {
  # Add "/" if it doesn't exist
  output_directory <- paste0(output_directory, "/")
}

# use standard function to find WT ORFs

#  run filter_custom_gtf, check if counts are present
if(!grepl("proteome_database_transcripts.gtf",gtf_path)){
  
  if (!is.null(tx_count_path)) {
    filter_custom_gtf(customgtf=gtf_path, organism=organism, tx_counts=tx_count_path, min_count=minimum_tx_count, outdir=output_directory)
    filtered_gtf <-paste0(output_directory, "proteome_database_transcripts.gtf") 
    
  } else {
    filter_custom_gtf(customgtf=gtf_path, organism=organism, outdir=output_directory)
    filtered_gtf <-paste0(output_directory, "proteome_database_transcripts.gtf")
    
  }
  
}else{
  
  filtered_gtf <-gtf_path
  
}

wt_orfome <- get_transcript_orfs(filteredgtf=filtered_gtf, organism=organism, orf_len=min_orf_length, find_UTR_5_orfs=find_5_orfs, find_UTR_3_orfs=find_3_orfs, referencegtf=reference_gtf,outdir=output_directory)


# generate custom genome
if (!is.null(ref_genome) && !is.null(vcf_file)) { # if a genome and vcf have been uploaded
  
  # define bash script command to inject variants into genome
  custom_genome_command <- paste0("bash bin/database_module/generate_custom_genome.sh -g  ", ref_genome, " -r ", reference_gtf, " -v ", vcf_file, " -o ", output_directory) #bin/database_module/generate_custom_genome.sh
  #print(custom_genome_command)
  system(custom_genome_command)
  
  # fetch variant protein sequences based on variants provided in VCF file
  genome_alt_hm_ht <- paste0(output_directory, basename(ref_genome) %>% str_replace(., ".fa", "_hm_ht.fa"))
  genome_alt_hm <- paste0(output_directory, basename(ref_genome) %>% str_replace(., ".fa", "_hm.fa"))
  # apply variant protein function
  get_variant_protein_seqs(wt_orfome=wt_orfome,custom_genome_hm =genome_alt_hm,custom_genome_hm_ht =genome_alt_hm_ht,  custom_gtf=paste0(output_directory, "proteome_database_transcripts.gtf"), organism=organism, outdir=output_directory, min_orf_length)
  
}







