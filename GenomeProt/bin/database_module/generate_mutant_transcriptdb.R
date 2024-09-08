# define functions
rm(list=ls())
#source("global.R")

library(optparse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(vcfR)
library(readr)



# option_list = list(
#   make_option(c("-r", "--reference"), type="character", default=NULL,
#               help="reference GTF", metavar="character"),
#   make_option(c("-v", "--variants"), type="character", default=NULL,
#               help="Varints in  tab delimited format", metavar="character"),
#   make_option(c("-t", "--transcript_db"), type="character", default=NULL,
#               help="Transcript database", metavar="character")
# )
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)


# reference_gtf <- opt$reference
# variants<- opt$variants
# transcript<-opt$transcript_db

genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
reference_gtf <- "/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/transcript.gtf"
vcf_file<-"/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/COSMIC.vcf"



#vcf file
vcf <- read.vcfR(vcf_file)
# Extract positions from the VCF file
chrom =vcf@fix[, "CHROM"]
ref_allele = vcf@fix[, "REF"]
alt_allele = vcf@fix[, "ALT"]
start = as.numeric(vcf@fix[, "POS"])
end = start+ (nchar(ref_allele) - 1)


variants_df<-data.frame(
  chrom = chrom,
  start = start,
  end = end,
  ref_allele = ref_allele,
  alt_allele = alt_allele
  
)



variants_df<-variants_df%>% unique()
write_tsv(variants_df,"/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/variants.tsv",col_names = FALSE)


#gtf
gtf <- import(reference_gtf)
exons <- gtf[gtf$type == "exon"]
exon_sequences <- getSeq(genome, exons)
mcols(exons)$sequence<-as.character(exon_sequences)


exon_with_sequences <- data.frame(
  seqnames = seqnames(exons),
  start = start(exons),
  end = end(exons),
  strand = strand(exons),
  transcript=mcols(exons)$transcript_id,
  transcript_type=mcols(exons)$transcript_type,
  exon=mcols(exons)$gene_id,
  gene_name=mcols(exons)$gene_name,
  sequence = mcols(exons)$sequence
)

write_tsv(exon_with_sequences,"/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/exons_seqs.tsv",col_names = FALSE)

#overlap using bedtools

command_bedtools<-system(paste0("bedtools intersect -loj -a ", "/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/variants.tsv"," -b ", "/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/exons_seqs.tsv"," | awk -F'\t' '{ if($6!=\".\") print $0 }' | uniq "),intern = TRUE)
  
df_bedtools <- read.table(text = command_bedtools, sep = "\t", header = FALSE) 

colnames(df_bedtools)<-c("chr","start","end","ref_allele","alt_allele","t_chr","e_start","e_end","strand","transcript_id","transcript_biotype","gene_id","gene_name","exon_seq")

df_bedtools<- df_bedtools%>% arrange(transcript_id)

compute_nucleotide_position <- function(start, e_start) {
  return(start - e_start+1)  # Position within the exon sequence
}

df_bedtools$nucleotide_position <- mapply(compute_nucleotide_position, df_bedtools$start, df_bedtools$e_start)


#reverse complement

reverse_complement <- function(seq,strand) {
  # Create a named vector for complement mapping
  complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  # Replace each nucleotide with its complement
  complemented_seq <- complement[unlist(strsplit(seq, split = ""))]
  
  # Reverse the sequence
  reverse_complement_seq <- paste(rev(complemented_seq), collapse = "")
  
  if (strand == "-") {
    reverse_complement_seq <- paste(rev(complemented_seq), collapse = "")
  } else if(strand == "+" ) {
    reverse_complement_seq <- paste(seq, collapse = "")
  }
  return(reverse_complement_seq)
}

df_bedtools$reverse_complement_exon_seq <- mapply(reverse_complement,df_bedtools$exon_seq, df_bedtools$strand)

extract_nucleotide <- function(seq, position) {
  # Convert 1-based position to 0-based for indexing
  position <- position
  if (position > nchar(seq) || position < 1) {
    return(NA)  # Return NA if the position is out of bounds
  }
  return(substr(seq, position, position))
}


df_bedtools$ref_nucleotide <- mapply(extract_nucleotide, df_bedtools$reverse_complement_exon_seq, df_bedtools$nucleotide_position)

df_bedtools<-df_bedtools%>%dplyr::select(chr,start,end,ref_allele,alt_allele,nucleotide_position,ref_nucleotide, everything())






