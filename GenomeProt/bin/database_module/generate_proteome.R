
# define functions

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
      mutate(total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
      dplyr::filter(total > as.numeric(min_count))
    
    # extract txnames
    tx_ids <- counts_filt$TXNAME
    
    # filter for these transcripts
    bambu_data <- bambu_data[mcols(bambu_data)$transcript_id %in% tx_ids]
    
  } 
  
  # remove scaffolds and weird chromosomes
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
  gene_query <- queryMany(unique(bambu_df$ensg_id), scopes="ensembl.gene", fields="symbol", species=as.character(organism),  returnall=TRUE)
  
  # make df
  gene_df <- as.data.frame(gene_query[["response"]])
  
  # if there was no name found, use original ID
  gene_df_formatted <- gene_df %>% 
    dplyr::mutate(gene_name = case_when(
      is.na(symbol) ~ query,
      !is.na(symbol) ~ symbol
    )) %>% 
    dplyr::select(query, gene_name)
  
  # merge results 
  bambu_merged <- merge(bambu_df, gene_df_formatted, by.x="ensg_id", by.y="query", all.x=T, all.y=F)
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

# fetch variant protein sequences based on VCF file
get_variant_protein_seqs <- function(custom_genome, custom_gtf, organism, wt_orfome, outdir, orf_len) {
  
  if (file.exists(custom_genome)) {
    # import filtered gtf as a txdb
    txdb <- makeTxDbFromGFF(custom_gtf)
    txs <- exonsBy(txdb, by=c("tx", "gene"), use.names=TRUE)
    # convert txdb to GRangesList
    txs_grl <- GRangesList(txs)
    alt_genome_fa <- FaFile(custom_genome)
    indexFa(custom_genome)
    
    # fetch variant transcript sequences
    mut_sequences <- GenomicFeatures::extractTranscriptSeqs(alt_genome_fa, txs_grl)
    
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
  
    # Fetch ORF nucleotide sequences
    orf_seqs <- GenomicFeatures::extractTranscriptSeqs(alt_genome_fa, ORFs)
    
    # convert nucleotide sequences to amino acid sequences
    orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,protein_sequence = orf_aa_seq, row.names=NULL)
    
    # remove special characters from protein sequence and add column with transcript id
    orf_aa_seq_df <- orf_aa_seq_df %>%
      mutate(transcript=str_replace(ORF_id, "_.*", "")) %>%
      mutate(protein_sequence=str_replace(protein_sequence,"\\*$", ""))
    
    # combine protein sequences with ORF genomic coordinates
    orf_aa_seq_df_genome_coord <- left_join(orf_aa_seq_df,orf_genome_coordinates,by=c("ORF_id"="names"))

    fasta_df_mut <- data.frame(
      transcript = names(mut_sequences),
      mut_seq = as.character(mut_sequences),
      stringsAsFactors = FALSE
    )
    
    #############
    
    # wild type
    
    # set organism
    if (organism == "human") {
      library(BSgenome.Hsapiens.UCSC.hg38)
      genomedb <- BSgenome.Hsapiens.UCSC.hg38
    } else if (organism == "mouse") {
      library(BSgenome.Mmusculus.UCSC.mm39)
      genomedb <- BSgenome.Mmusculus.UCSC.mm39
    } else if (organism == "celegans") {
      library(BSgenome.Celegans.UCSC.ce11)
      genomedb <- BSgenome.Celegans.UCSC.ce11
    } else if (organism == "drosophila") {
      library(BSgenome.Dmelanogaster.UCSC.dm6)
      genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
    } else if (organism == "rat") {
      library(BSgenome.Rnorvegicus.UCSC.rn7)
      genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
    } else if (organism == "zebrafish") {
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
  
    # function to calculate the number of variable nucleiotides
    count_variable_nucleotides <- function(mut_seq, wt_seq) {
      sum(strsplit(mut_seq, "")[[1]] != strsplit(wt_seq, "")[[1]])
    }
    
    # apply the function to the dataframe
    fasta_df_merged$variable_nucleotides <- mapply(count_variable_nucleotides, fasta_df_merged$mut_seq, fasta_df_merged$wt_seq)
    
    # remove sequences with no variable nucleotides
    fasta_df_merged_mutant <- fasta_df_merged %>%
      filter(variable_nucleotides != 0) %>% 
      dplyr::select(transcript,mut_seq) %>%
      unique()
    
    # subset ORFs for variant transcripts
    orf_aa_seq_df_genome_coord_filtered <- orf_aa_seq_df_genome_coord %>%
      filter(transcript %in% fasta_df_merged_mutant$transcript)
    
    # adding the transcript_id column by splitting ORF_id and taking the first part
    orfome_wt_df <- wt_orfome %>% 
      mutate(transcript = sapply(strsplit(ORF_id, "_"), `[`, 1)) %>%
      dplyr::select(ORF_sequence,transcript)
    
    # merge wild type protein sequences from ORFome_aa dataframe and dataframe (orf_aa_seq_df_genome_coord_filtered) containing variants and wild type sequences based on transcript ids
    mutated_transcript_orfs_combined <- inner_join(orf_aa_seq_df_genome_coord_filtered,orfome_wt_df,by = "transcript") %>% 
      mutate(ORF_sequence=str_replace(ORF_sequence,"\\*", ""))
    
    # remove the wild type sequences from protein_sequence column of orf_aa_seq_df_genome_coord_filtered
    mutated_transcript_orfs_combined_flt <- mutated_transcript_orfs_combined %>%
      group_by(transcript) %>%
      filter(protein_sequence!=ORF_sequence) %>%
      ungroup()
    
    # function to perform global alignment and calculate identical residues
    calculate_identity_info <- function(orf_seq, protein_seq) {
      # perform global alignment
      alignment <- pairwiseAlignment(orf_seq, protein_seq, substitutionMatrix = "BLOSUM62", type = "global")
      
      # extract the number of identical residues
      identical_residues <- nmatch(alignment)
      
      # get the total alignment length (length of the aligned sequences)
      alignment_length <- nchar(pattern(alignment))
      
      # calculate percent identity
      percent_identity <- (identical_residues / alignment_length) * 100
      
      return(c(identical_residues,percent_identity))
    }
    
    mutated_transcript_orfs_combined_flt_identity <- mutated_transcript_orfs_combined_flt %>%
      rowwise() %>%
      mutate(identity_info = list(calculate_identity_info(ORF_sequence, protein_sequence))) %>%
      mutate(identical_residues = identity_info[[1]], percent_identity = identity_info[[2]]) %>%
      dplyr::select(-identity_info) %>%
      ungroup()

    mutated_transcript_orfs_combined_flt_90per_identity <- mutated_transcript_orfs_combined_flt_identity %>%
      dplyr::filter(percent_identity >= 90) %>%
      dplyr::filter(percent_identity < 100)
    
    variant_protein_seqs <- mutated_transcript_orfs_combined_flt_90per_identity %>%
      mutate(orf_coodinates=paste0(chr,":",start,"-",end)) %>%
      dplyr::select(transcript,protein_sequence,orf_coodinates) %>%
      unique()
    
    print(mutated_transcript_orfs_combined_flt_90per_identity)
    
    # export protein seqs for python script
    write_tsv(variant_protein_seqs, paste0(outdir, "/Mutant_ORFome_aa.txt"))

  } else {
    
    print("File does not exist")
  
  }
 
}

# generate FASTA of transcript sequences
get_transcript_orfs <- function(filteredgtf, organism, orf_len=30, find_UTR_5_orfs=FALSE, find_UTR_3_orfs=FALSE, referencegtf, outdir) {
  
  # set organism
  if (organism == "human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
  } else if (organism == "mouse") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39
  } else if (organism == "celegans") {
    library(BSgenome.Celegans.UCSC.ce11)
    genomedb <- BSgenome.Celegans.UCSC.ce11
  } else if (organism == "drosophila") {
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
  } else if (organism == "rat") {
    library(BSgenome.Rnorvegicus.UCSC.rn7)
    genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
  } else if (organism == "zebrafish") {
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
  
  # ORFik function to find ORFs in GRangesList
  apply_orfik <- function(grl, orfik_min_length, orfik_max_length, orfik_type) {
    
    # extract transcript nt sequence
    tx_seqs <- extractTranscriptSeqs(genomedb, grl)
    
    # ORFik
    ORFs <- findMapORFs(grl,
                        tx_seqs, 
                        groupByTx = FALSE,
                        longestORF = orfik_type, 
                        minimumLength = as.numeric(orfik_min_length), 
                        startCodon = "ATG",
                        stopCodon = stopDefinition(1))
    
    # unlist GRL
    ORFs_unlisted <- unlist(ORFs) %>% as_tibble()
    
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
    
    # convert the nucleotide sequences to amino acid sequences
    orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,ORF_sequence = orf_aa_seq, row.names=NULL) 
    
    # merge with coordinates
    orf_aa_seq_df_genomic_coordinates <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = "ORF_id")
    
    return(orf_aa_seq_df_genomic_coordinates)
    
  }
  
  # ORF discovery
  
  # set ORF max length to large number (to disable)
  # set longestORF to FALSE to ensure we identify CDS
  all_ORFs <- apply_orfik(txs_grl, as.numeric(orf_len), 100000, FALSE)
  
  # create tmp copy
  tmp <- all_ORFs
  
  # extract 5UTR ORFs
  if (find_UTR_5_orfs == TRUE) {
    
    # extract 5' UTR regions from ref gtf
    utrs5 <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs5_filtered <- utrs5[names(utrs5) %in% names(txs)]
    
    # ORF max length is now the main ORF min length
    # only keep longest UTR ORFs
    utr5_ORFs <- apply_orfik(utrs5_filtered, 10, as.numeric(orf_len), TRUE)
    
    # update tmp
    tmp <- rbind(tmp, utr5_ORFs)
    
  }
  
  # extract 3UTR ORFs
  if (find_UTR_3_orfs == TRUE) {
    
    # extract 3' UTR regions from ref gtf
    utrs3 <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs3_filtered <- utrs3[names(utrs3) %in% names(txs)]
    
    # ORF max length is now the main ORF min length
    # only keep longest UTR ORFs
    utr3_ORFs <- apply_orfik(utrs3_filtered, 10, as.numeric(orf_len), TRUE)
    
    # update tmp
    tmp <- rbind(tmp, utr3_ORFs)
    
  }
  
  # rename final ORFs with new numerical IDs
  combined <- tmp %>% 
    separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
    group_by(transcript_id) %>% 
    mutate(tx_id_number = row_number()) %>% 
    ungroup()
  
  # apply new names
  combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
  combined$transcript_id <- NULL
  combined$tx_id_number <- NULL
  
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
              help="reference user GTF", metavar="character"),
  make_option(c("-G", "--genome"), type="character", default=NULL,
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

#  run filter_custom_gtf, check if counts are present
if (!is.null(tx_count_path)) { # if counts have been uploaded
 filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, organism=organism, tx_counts=tx_count_path, min_count=minimum_tx_count, outdir=output_directory)
} else { # no counts
 filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, organism=organism, outdir=output_directory)
}

# use standard function to find WT ORFs
wt_orfome <- get_transcript_orfs(filteredgtf=paste0(output_directory, "/proteome_database_transcripts.gtf"), organism=organism, orf_len=min_orf_length, find_UTR_5_orfs=find_5_orfs, find_UTR_3_orfs=find_3_orfs, referencegtf=reference_gtf, outdir=output_directory)

# generate custom genome
if (!is.null(ref_genome) && !is.null(vcf_file)) { # if a genome and vcf have been uploaded
  
  # define bash script command to inject variants into genome
  # supply conda env here if local installation
  custom_genome_command <- paste0("bash bin/database_module/generate_custom_genome.sh -g  ", ref_genome, " -r ", reference_gtf, " -v ", vcf_file, " -o ", output_directory, "/")
  system(custom_genome_command)
 
  # fetch variant protein sequences based on variants provided in VCF file
  genome_alt <- paste0(output_directory, "/", basename(ref_genome) %>% str_replace(., ".fa", "_alt.fa"))
 
  # apply variant protein function
  get_variant_protein_seqs(custom_genome=genome_alt, custom_gtf=paste0(output_directory, "/proteome_database_transcripts.gtf"), organism=organism, wt_orfome=wt_orfome, outdir=output_directory, min_orf_length)

}
 
 
 
 

 


