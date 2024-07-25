
run_bambu_function <- function(bam_file_list, gtf, genome) { 

  system(paste0("mkdir bambu_output"))
  
	bambuAnnotations <- prepareAnnotations(gtf)
	se <- bambu(reads = bam_file_list, annotations = bambuAnnotations, genome = genome)
	writeBambuOutput(se, path = "bambu_output")
	
	tx_data <- as.data.frame(mcols(se))
	tx_data <- as.data.frame(apply(tx_data, 2, as.character))
	tx_data <- tx_data %>% dplyr::filter(novelTranscript == "TRUE")
	
	write.csv(tx_data, "bambu_output/novel_transcript_classes.csv", row.names=F, quote=F)
	
	system(paste0("mv bambu_output/counts_transcript.txt bambu_output/bambu_transcript_counts.txt"))
	system(paste0("mv bambu_output/extended_annotations.gtf bambu_output/bambu_transcript_annotations.gtf"))
	
}
