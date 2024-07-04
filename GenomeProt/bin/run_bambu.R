library("optparse")

run_bambu_function <- function(bam_dir, gtf, genome) { 

  bam_files <- as.vector(list.files(bam_dir, pattern = "\\.bam", full.names=TRUE))
  
	bambuAnnotations <- prepareAnnotations(gtf)
	se <- bambu(reads = bam_files, annotations = bambuAnnotations, genome = genome)
	tx_data <- as.data.frame(mcols(se))
	tx_data <- apply(tx_data, 2, as.character)
	
	write.table(tx_data, "bambu_tx_classes.txt", row.names=F, quote=F, sep="\t")
	writeBambuOutput(se, path = '.')
	
}

option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="BAM file directory", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Reference annotation file (GTF/GFF)", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Genome fasta file", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$bam)){
  
  print_help(opt_parser)
  stop("Incorrect arguments", call.=FALSE)

} else if (file.exists(opt$gtf) & file.exists(opt$fasta) & dir.exists(opt$bam)) {
  
  library(bambu)
  library(dplyr)

  run_bambu_function(opt$bam,opt$gtf,opt$fasta)
  
} else { 
  
  
  print_help(opt_parser)
  stop("Incorrect arguments", call.=FALSE)

}
