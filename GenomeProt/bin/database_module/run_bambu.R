library(optparse)
library(bambu)
library(stringr)
library(dplyr)
library(Rsamtools)


run_bambu_function <- function(bam_file_list, gtf, organism, output_directory) { 
  
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
  
  bambuAnnotations <- prepareAnnotations(gtf)
  se <- bambu(reads = bam_file_list, annotations = bambuAnnotations, genome = genomedb, verbose = TRUE)
  writeBambuOutput(se, path = output_directory)
  
  tx_data <- as.data.frame(mcols(se))
  tx_data <- as.data.frame(apply(tx_data, 2, as.character))
  tx_data <- tx_data %>% dplyr::filter(novelTranscript == "TRUE")
  
  write.csv(tx_data, paste0(output_directory, "/novel_transcript_classes.csv"), row.names=F, quote=F)
  
}



option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="BAM file directory (optional)", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Transcript annotation file (GTF/GFF)", metavar="character"),
  make_option(c("-s", "--species: mouse,celegans,drosophila,rat,zebrafish"), type="character", default=NULL,
              help="Organism: ", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


    
  if (!is.null(opt$bam) && file.exists(opt$bam) && 
      !is.null(opt$gtf) && file.exists(opt$gtf) && 
      !is.null(opt$species)){
      
          # create list of BAMs
          bam_files <- list.files(path = opt$bam, "\\.bam$", full.names = TRUE)
          bam_file_list <- Rsamtools::BamFileList(bam_files)
          print(bam_file_list)
          
          # remove bam extension
          bam_files_names <- list.files(path =opt$bam, "\\.bam$", full.names = FALSE)
          bam_file_names <- str_remove(bam_files_names, ".bam")
          # rename list to original names
          names(bam_file_list) <- bam_file_names
          print(bam_file_list)
          
        }


# run bambu function from R/ dir
run_bambu_function(bam_file_list, opt$gtf, opt$species, opt$outdir)

