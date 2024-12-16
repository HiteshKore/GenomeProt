library("optparse")

Bambu<- function (bam,gtf,genome.fa,outdir) { 
	bambuAnnotations <- prepareAnnotations(gtf)
	se <-bambu(reads=bam, annotations = bambuAnnotations, genome = genome.fa)
	tx_data <- as.data.frame(mcols(se))
	tx_data <- apply(tx_data,2,as.character)
	system(paste0('mkdir ',outdir))
	write.table(tx_data, paste0(outdir, "/bambu_tx_classes.txt"), row.names=F, quote=F, sep="\t")
	writeBambuOutput(se, path = outdir)

}


option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="BAM file directory", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Transcript annotation file (GTF/GFF)", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Genome fasta file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Incorrect arguments",call.=FALSE)
  
  
  
}else if ( file.exists(opt$gtf) & file.exists(opt$fasta) & dir.exists(opt$bam)){
  
  
  library("bambu")
  library(stringr)
  library(dplyr)

  files=list.files(opt$bam)
  bam_files=str_subset(files,".bam")
  Bam_files_path=paste0(opt$bam,"/",bam_files)

	if(length(Bam_files_path)==1){
	
		Bambu(Bam_files_path,opt$gtf,opt$fasta,opt$outdir)       
	 } else if (length(Bam_files_path)>1){

    Bambu(Bam_files_path,opt$gtf,opt$fasta,opt$outdir)

	}



}else{
 print_help(opt_parser)
 stop("Incorrect arguments",call.=FALSE)

}
