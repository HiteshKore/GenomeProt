
run_bambu_function <- function(bam_file_list, gtf, organism, threads=1) { 

  system(paste0("mkdir bambu_output"))
  
  if (organism=="mouse") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39
    
  } else if (organism=="human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
  }
  
	bambuAnnotations <- prepareAnnotations(gtf)
	se <- bambu(reads = bam_file_list, annotations = bambuAnnotations, genome = genomedb)
	writeBambuOutput(se, path = "bambu_output")
	
	tx_data <- as.data.frame(mcols(se))
	tx_data <- as.data.frame(apply(tx_data, 2, as.character))
	tx_data <- tx_data %>% dplyr::filter(novelTranscript == "TRUE")
	
	write.csv(tx_data, "bambu_output/novel_transcript_classes.csv", row.names=F, quote=F)
	
	system(paste0("mv bambu_output/counts_transcript.txt bambu_output/bambu_transcript_counts.txt"))
	system(paste0("mv bambu_output/extended_annotations.gtf bambu_output/bambu_transcript_annotations.gtf"))
	
}
