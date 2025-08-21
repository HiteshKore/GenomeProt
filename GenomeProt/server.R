
library(shiny)
library(shinyjs)
library(reticulate)

conda_path=conda_binary() #sometimes r environment doesn't recognize the conda path

# internal server functions
fastq_server <- function(input, output, session) {
  
  # to start with fastqs a user must upload a genome, reference gtf, fastq files
  req(input$user_reference_genome$datapath, input$user_reference_gtf$datapath, input$user_fastq_files$datapath)  # required
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_bam <- paste0(session_id, "/mapping_output")
  # create output dir
  system(paste0("mkdir ", outdir_bam))
  
  # optionally export genome from R
  # if (organism == "human") {
  #   library(BSgenome.Hsapiens.UCSC.hg38)
  #   genomedb <- BSgenome.Hsapiens.UCSC.hg38
  # } else if (organism == "mouse") {
  #   library(BSgenome.Mmusculus.UCSC.mm39)
  #   genomedb <- BSgenome.Mmusculus.UCSC.mm39
  # } else if (organism == "celegans") {
  #   library(BSgenome.Celegans.UCSC.ce11)
  #   genomedb <- BSgenome.Celegans.UCSC.ce11
  # } else if (organism == "drosophila") {
  #   library(BSgenome.Dmelanogaster.UCSC.dm6)
  #   genomedb <- BSgenome.Dmelanogaster.UCSC.dm6
  # } else if (organism == "rat") {
  #   library(BSgenome.Rnorvegicus.UCSC.rn7)
  #   genomedb <- BSgenome.Rnorvegicus.UCSC.rn7
  # } else if (organism == "zebrafish") {
  #   library(BSgenome.Drerio.UCSC.danRer11)
  #   genomedb <- BSgenome.Drerio.UCSC.danRer11
  # }
  # export(genomedb, "genome.fasta.gz", verbose=T, compress=T, format="fasta")
  
  # create dataframe with sample details and add prefix column 
  user_fastq_files_df <- input$user_fastq_files %>% 
    mutate(file_prefix = str_replace_all(name, c("\\.fastq\\.gz$" = "", "\\.fastq$" = "", "\\.fq$" = "", "\\.fa$" = "", "\\.fasta$" = "")))
  
  print(user_fastq_files_df)
  
  # map reads
  if (input$sequencing_type == "long-read") { # for long reads, use minimap2
    
    # define index file name for genome
    index_file <- paste0(outdir_bam, "/", str_replace_all(input$user_reference_genome$name, c("\\.fa$" = "", "\\.fasta$" = "")), ".mmi")
    
    # define minimap2 command to index genome
    minimap2_index_command <- paste0(conda_path," run -n GenomeProt_env minimap2 -ax splice:hq -d ", index_file, " ", input$user_reference_genome$datapath)
    
    print(minimap2_index_command)
    
    # run minimap2 indexing
    system(minimap2_index_command)
    
    # for each fastq file
    for (i in 1:nrow(user_fastq_files_df)) {
      
      fastq_file <- user_fastq_files_df$datapath[i] # set file path
      file_prefix <- user_fastq_files_df$file_prefix[i] # get file prefix
      
      # define command to map fastq file reads to the indexed genome file with minimap2
      minimap2_command <- paste0(conda_path," run -n GenomeProt_env minimap2 -t ", input$user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", index_file, " ", fastq_file, " | samtools view -bh | samtools sort -@ ", input$user_threads, " -o ", outdir_bam, "/", file_prefix, ".bam")
      
      print(minimap2_command)
      
      # run mapping 
      system(minimap2_command)
      
    }
    
  } else if (input$sequencing_type == "short-read") { # for short reads, use salmon
    
    # generate salmon index
    if ( grepl("\\.gz$", input$user_reference_genome$datapath) & grepl("\\.gz$", input$transcriptome_file$datapath)) {
      
      command_generate_decoy <- paste0("bash -c \"grep '^>' <(gunzip -c ", input$user_reference_genome$datapath, ") | cut -d ' ' -f 1 >", outdir_bam, "/decoys.txt",'\"')
      command_sed <- paste0("sed -i -e 's/>//g' ", outdir_bam, "/decoys.txt")
      command_ref_file <- paste0("cat ", input$transcriptome_file$datapath, " ", input$user_reference_genome$datapath, " >", outdir_bam, "/gentrome.fa.gz")
      command_index <- paste0(conda_path," run -n GenomeProt_env salmon index -t ", outdir_bam, "/gentrome.fa.gz -d ", outdir_bam, "/decoys.txt -p ", input$user_threads, " -i ", outdir_bam, "/salmon_index --gencode")
      
    } else {
      
      command_generate_decoy <- paste0("grep '^>' ", input$user_reference_genome$datapath, " | cut -d \" \" -f 1 >", outdir_bam, "/decoys.txt")
      command_sed <- paste0("sed -i -e 's/>//g' ", outdir_bam, "/decoys.txt")
      command_ref_file <- paste0("cat ", input$transcriptome_file$datapath, " ", input$user_reference_genome$datapath, " >", outdir_bam, "/gentrome.fa")
      command_index <- paste0(conda_path," run -n GenomeProt_env salmon index -t ", outdir_bam, "/gentrome.fa -d ", outdir_bam, "/decoys.txt -p ", input$user_threads, " -i ", outdir_bam, "/salmon_index --gencode")
      
    }
    
    # generate index
    system(command_generate_decoy)
    system(command_sed)
    system(command_ref_file)
    system(command_index)
    
    # determine paired end or single end
    paired_end <- list(list(R1 = NULL, R2 = NULL))
    single_end <- list()
    
    # for each fastq file
    for (i in 1:nrow(user_fastq_files_df)) {
      
      file_name <- user_fastq_files_df$name[i]
      
      if (grepl("_R1", file_name)) {
        
        base_name <- sub("_R1(\\.fastq(\\.gz)?)$", "", file_name)
        paired_end[[base_name]]$R1 <- user_fastq_files_df$datapath[i]
        
      } else if (grepl("_R2", file_name)) {
        
        base_name <- sub("_R2(\\.fastq(\\.gz)?)$", "", file_name)
        paired_end[[base_name]]$R2 <- user_fastq_files_df$datapath[i]
        
      }
      
      if(!grepl("_R1", file_name) && !grepl("_R2", file_name)) {
        
        base_name <- sub("(\\.fastq(\\.gz)?)$", "", file_name)
        single_end[[base_name]] <- user_fastq_files_df$datapath[i]
        
      } 
    }
    
    # if there are paired end reads
    if (length(paired_end) > 0) {
      
      for (base_name in names(paired_end)) {
        # access R1 and R2 elements
        R1_path <- paired_end[[base_name]]$R1
        R2_path <- paired_end[[base_name]]$R2
        
        if (!is.null(R1_path) && !is.null(R2_path)) {
          
          # perform operations with R1 and R2
          cat("Base Name:", base_name, "\t","R1 Path:", R1_path, "R2 Path:", R2_path, "\n")
          
          # define salmon command
          command_salmon <- paste0(conda_path," run -n GenomeProt_env salmon quant -i ", outdir_bam,"/salmon_index -p ", input$user_threads ," -l A -1 ", R1_path," -2 ", R2_path, " --validateMappings -o ", outdir_bam,"/", base_name)
          print(command_salmon)
          # run salmon
          system(command_salmon)
          
        }
      }
    }
    
    # if there are single end reads
    if (length(single_end) > 0) {
      for (base_name in names(single_end)) {
        
        print(single_end[[base_name]])
        
        # define salmon command
        command_salmon <- paste0(conda_path," run -n GenomeProt_env salmon quant -i ", outdir_bam,"/salmon_index -p ", input$user_threads ," -l A -r ", single_end[[base_name]]," --validateMappings -o ", outdir_bam,"/", base_name)
        print(command_salmon)
        # run salmon
        system(command_salmon)
        
      }
    }
    
    # create count matrix
    samples <- c(names(single_end),names(paired_end))
    samples <- samples[samples != ""]
    print(samples)
    # set path to salmon quant files
    files <- file.path(outdir_bam, samples, "quant.sf")
    print(files)
    names(files) <- samples
    
    # check if the files exist
    print(all(file.exists(files)))
    
    # use tximport to import salmon quantification files
    txi <- tximport(files, type = "salmon", txOut = TRUE)
    print(input$user_reference_gtf$datapath)
    # import gtf to add gene information
    gtf_data <- rtracklayer::import(input$user_reference_gtf$datapath, format = "gtf")
    
    # convert GTF data to a data frame
    gtf_df <- as.data.frame(gtf_data)
    
    # filter for relevant columns
    transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
    colnames(transcript_gene_info) <- c("TXNAME","GENEID")
    
    # convert counts object to df
    count_df <- as.data.frame(txi$counts)
    # set col names
    count_df <-  count_df %>% mutate(TXNAME = rownames(count_df)) %>% dplyr::select(TXNAME, everything())
    rownames(count_df) <- NULL
    # merge tx counts and gene info
    count_df_merged <- left_join(count_df, transcript_gene_info, by="TXNAME")
    
    count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, everything())
    
    # export short-read count data
    write_tsv(count_df_merged, file = paste0(outdir_bam, "/transcript_counts.txt"), escape = "none", col_names = TRUE)
    
  }
}

bam_server <- function(input, output, session) {
  
  req(input$user_reference_genome_bam$datapath,input$user_reference_gtf$datapath )  # BAMs and ref GTF required
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_bam <- paste0(session_id, "/mapping_output")
  # create outpir dir
  system(paste0("mkdir ", outdir_bam))
  
  # generate reference.fa
  # check if file is compressed
  if (grepl("\\.gz$", input$user_reference_genome_bam$datapath)) {
    
    command_decompress <- paste0(conda_path," run -n GenomeProt_env gzip -c -d ", input$user_reference_genome_bam$datapath, " >", outdir_bam, "/genome.fa")
    print(command_decompress)
    system(command_decompress)
    command_gffread <- paste0(conda_path," run -n GenomeProt_env gffread -w ", outdir_bam, "/transcript.fa -g ", outdir_bam, "/genome.fa ", input$user_reference_gtf$datapath)
    
  } else {
    
    command_gffread <- paste0(conda_path," run -n GenomeProt_env gffread -w ", outdir_bam, "/transcript.fa -g ", input$user_reference_genome_bam$datapath, " ", input$user_reference_gtf$datapath)
    
  }
  
  print(command_gffread)
  system(command_gffread) 
  
  # create df of bam file names
  user_bam_files_df <- input$user_bam_files %>%
    mutate(file_prefix = sub("\\.bam$", "", name))
  
  print(user_bam_files_df)
  
  # for each bam file
  for (i in 1:nrow(user_bam_files_df)) {
    
    bam_file <- user_bam_files_df$datapath[i]
    file_prefix <- user_bam_files_df$file_prefix[i]
    
    # define salmon command
    command_salmon <- paste0(conda_path," run -n GenomeProt_env salmon quant -t ", outdir_bam, "/transcript.fa -p ", input$user_threads , " -l A  -a ", bam_file , " -o ", outdir_bam, "/", file_prefix)
    
    print(command_salmon)
    system(command_salmon)
    
  }
  
  # create count matrix
  samples <- user_bam_files_df$file_prefix
  # set path to salmon quant files
  files <- file.path(outdir_bam, samples, "quant.sf")
  names(files) <- samples
  
  # check if the files exist
  print(all(file.exists(files)))
  
  # use tximport to import salmon quantification files
  txi <- tximport(files, type = "salmon", txOut = TRUE)
  
  # import gtf to add gene information
  gtf_data <- rtracklayer::import(input$user_reference_gtf$datapath, format = "gtf")
  
  # convert GTF data to a data frame
  gtf_df <- as.data.frame(gtf_data)
  
  # filter for relevant columns
  transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
  colnames(transcript_gene_info) <- c("TXNAME","GENEID")
  
  # convert counts object to df
  count_df <- as.data.frame(txi$counts)
  # set col names
  count_df <-  count_df %>% mutate(TXNAME = rownames(count_df)) %>% dplyr::select(TXNAME, everything())
  rownames(count_df) <- NULL
  # merge tx counts and gene info
  count_df_merged <- left_join(count_df, transcript_gene_info, by="TXNAME")
  
  count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, everything())
  
  # export short-read count data
  write_tsv(count_df_merged, file = paste0(outdir_bam, "/transcript_counts.txt"), escape = "none", col_names = TRUE)
  
}

bambu_server <- function(input, output, session) {
  
  # store session ID
  session_id <- session$token
  # set BAM dir
  outdir_bam <- paste0(session_id, "/mapping_output")
  # set output dir
  outdir_bambu <- paste0(session_id, "/bambu_output")
  # create output dir
  system(paste0("mkdir ", outdir_bambu))
  
  
   if (input$input_type == "bam_input") { # if user input bam files
     req(input$user_bam_files$datapath, input$user_reference_gtf$datapath, input$organism)  # required
     
     #print(input$user_bam_files)
     
     bam_df=as.data.frame(input$user_bam_files)
     
     #print(bam_df)
     
     bamdir=dirname(input$user_bam_files$datapath)
     
     # Specify new filenames (you can adjust based on logic)
     new_names <- paste0(bamdir[1],"/", input$user_bam_files$name)
     print(new_names)
     # Rename the files
     file.rename(input$user_bam_files$datapath, new_names)
     print(input$user_bam_files)
    
     command_bambu= paste0(conda_path," run -n GenomeProt_env && Rscript bin/database_module/run_bambu.R -b ", bamdir[1],
     "/ -g ",input$user_reference_gtf$datapath," -o ", outdir_bambu," -s ", input$organism)
     
     

   
    } else if (input$input_type == "fastq_input") { # if user input fastq files
  
     command_bambu= paste0(conda_path, " run -n GenomeProt_env && Rscript bin/database_module/run_bambu.R -b ", outdir_bam,
                           " -g ",input$user_reference_gtf$datapath," -o ", outdir_bambu," -s ", input$organism)
     
    
  }
  
  #execute bambu command
  print(command_bambu)
  system(command_bambu)
  
  # rename bambu output files
  system(paste0("mv ",  outdir_bambu, "/counts_transcript.txt ", outdir_bambu, "/transcript_counts.txt"))
  system(paste0("mv ",  outdir_bambu, "/extended_annotations.gtf ", outdir_bambu, "/bambu_transcript_annotations.gtf"))
  
  # define gffcompare command for bambu gtf
  #command_gff_compare <- paste0("source activate IsoLamp; gffcompare -r ", input$user_reference_gtf$datapath, " ", outdir_bambu, "/bambu_transcript_annotations.gtf")
  command_gff_compare <- paste0(conda_path," run -n GenomeProt_env gffcompare -r ", input$user_reference_gtf$datapath, " ", outdir_bambu, "/bambu_transcript_annotations.gtf")
  print(command_gff_compare)
  # run gffcompare
  system(command_gff_compare)
  
  # rename gffcompare output files
  system(paste0("mv ", outdir_bambu, "/gffcmp.bambu_transcript_annotations.gtf.tmap ", outdir_bambu, "/gffcompare.tmap.txt"))
  system(paste0("rm gffcmp*"))
  
}

database_server <- function(input, output, session) {
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_db <- paste0(session_id, "/database_output")
  # create output dir
  system(paste0("mkdir ", outdir_db))
  
  # set input file type
  if (input$input_type == "gtf_input") {
    req(input$user_gtf_file, input$user_reference_gtf) # GTFs required
    db_gtf_file <- input$user_gtf_file$datapath
    db_counts_file <- input$user_tx_count_file$datapath
  } else if (input$input_type %in% c("bam_input", "fastq_input")) {
    if (input$sequencing_type == "long-read") {
      db_gtf_file <- paste0(session_id, "/bambu_output/bambu_transcript_annotations.gtf")
      db_counts_file <- paste0(session_id, "/bambu_output/transcript_counts.txt")
    } else if (input$sequencing_type == "short-read") {
      db_gtf_file <- input$user_reference_gtf$datapath
      db_counts_file <- paste0(session_id, "/mapping_output/transcript_counts.txt")
    }
  }
  
  # conditionally add counts file argument
  counts_arg <- if (!is.null(db_counts_file) && nzchar(db_counts_file)) {
    paste0(" -c ", db_counts_file)
  } else {
    ""
  }
  
  # conditionally add genome input argument
  genome_arg <- if (!is.null(input$user_reference_genome$datapath)) {
    paste0(" -G ", input$user_reference_genome$datapath)
  } else {
    ""
  }
  
  # construct the command
  command_generate_proteome <- paste0(
    conda_path," run -n GenomeProt_env && Rscript bin/database_module/generate_proteome.R",
    " -g ", db_gtf_file,
    " -r ", input$user_reference_gtf$datapath,
    counts_arg, # include counts file only if provided
    " -m ", input$minimum_tx_count,
    " -o ", input$organism,
    " -l ", input$min_orf_length,
    " -u ", input$user_find_utr_5_orfs,
    " -d ", input$user_find_utr_3_orfs,
    if (!is.null(input$user_vcf_file)) paste0(" -v ", input$user_vcf_file$datapath) else "",
    genome_arg, # include genome file only if provided
    " -s ", outdir_db
  )
  
  # print command
  print(command_generate_proteome)
  
  # run command
  system(command_generate_proteome)
  
  print("Generated ORFs")
  
  # set reference protein database per organism 
  if (input$organism == "HUMAN") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "MOUSE") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  } else if (input$organism == "CAEEL") {
    ref_proteome <- "data/openprot_uniprotDb_c_elegans.txt"
  } else if (input$organism == "DROME") {
    ref_proteome <- "data/openprot_uniprotDb_drosophila.txt"
  } else if (input$organism == "RAT") {
    ref_proteome <- "data/openprot_uniprotDb_rat.txt"
  } else if (input$organism == "DANRE") {
    ref_proteome <- "data/openprot_uniprotDb_zebrafish.txt"
  }
  
  # run python script using conda env
   
  if (!is.null(input$user_vcf_file)) { # if there is a VCF file uploaded
    command_annotate_proteome <- paste0(conda_path," run -n GenomeProt_env --no-capture-output python bin/database_module/annotate_proteome.py ", input$user_reference_gtf$datapath, " ", ref_proteome, " ", outdir_db, "/ORFome_aa.txt ", outdir_db, "/proteome_database_transcripts.gtf ", outdir_db, " ", input$database_type, " ", input$min_orf_length, " ", paste0(outdir_db, "/Mutant_ORFome_aa.txt "),input$organism )
    print(command_annotate_proteome)
  } else { # if no VCF file uploaded
    command_annotate_proteome <- paste0(conda_path," run -n GenomeProt_env --no-capture-output python bin/database_module/annotate_proteome.py ", input$user_reference_gtf$datapath, " ", ref_proteome, " ", outdir_db, "/ORFome_aa.txt ", outdir_db, "/proteome_database_transcripts.gtf ", outdir_db, " ", input$database_type, " ", input$min_orf_length, " None ", input$organism)
    print(command_annotate_proteome)
  }
  
  # run python script to create proteome fasta
  system(command_annotate_proteome)
  
  print("Annotated proteome")
  
  # get top level directory
  top_level_dir <- getwd()
  
  # zip all results files depending on input types
  if (file.exists(paste0(outdir_db, "/proteome_database.fasta")) && file.exists(paste0(outdir_db, "/proteome_database_transcripts.gtf")) && !file.exists(paste0(outdir_db, "/orf_temp.txt"))) {
    if (input$input_type == "fastq_input" & input$sequencing_type == "long-read") {
      bam_files <- list.files(path = paste0(session_id, "/mapping_output"), "\\.bam$", full.names = TRUE)
      files_to_zip_db <- c(bam_files, "../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    } else if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      files_to_zip_db <- c("../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    } else if (input$sequencing_type == "short-read") {
      files_to_zip_db <- c("../mapping_output/transcript_counts.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    } else if (input$input_type == "gtf_input") {
      files_to_zip_db <- c("proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    }
    
    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_db <- file.path("../database_results.zip")
    
    # temp change the working dir to outdir_db
    tmp_wd <- setwd(outdir_db)
    
    # zip files
    zip(zipfile = zipfile_path_db, files = files_to_zip_db)
    
    # change back to starting wd
    setwd(top_level_dir)
    
  }
  
}

proteomics_server <- function(input, output, session) {
  
  # req(input$user_mm_data, input$user_mm_fasta)
  # 
  # # get directory path
  # dir_path <- dirname(input$user_mm_data$datapath[1])
  # print(dir_path)
  # 
  # mass_spec_names <- c(input$user_mm_data$name)
  # print(mass_spec_names)
  # 
  # # new file paths
  # new_file_paths <- file.path(dir_path, mass_spec_names)
  # 
  # # rename files
  # file.rename(input$user_mm_data$datapath, new_file_paths)
  # 
  # # copy config files
  # # replace all lines:
  # # 'MaxThreadsToUsePerFile = 3'
  # # with user set threads
  # 
  # #system(paste0("source activate mm_env; metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  # system(paste0("metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  # 
  # # check files exist
  # if (file.exists("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv") && file.exists("proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")) {
  #   # create a zip file with results
  #   files_to_zip <- c("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv", "proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
  #   zipfile_path <- "proteomics_output/proteomics_results.zip"
  #   zip(zipfile = zipfile_path, files = files_to_zip)
  # }
}

integration_server <- function(input, output, session) {
  
  req(input$user_proteomics_file, input$user_post_gtf_file, input$user_metadata_file)  # GTF is required
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_integ <- paste0(session_id, "/integ_output")
  # create output dir
  system(paste0("mkdir ", outdir_integ))
  
  # run Rscript
  system(paste0( conda_path," run -n GenomeProt_env && Rscript bin/integration_module/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -m ", input$user_metadata_file$datapath, " -g ", input$user_post_gtf_file$datapath, " -s ", outdir_integ))
  
  # get the top level dir
  top_level_dir <- getwd()
  
  system(paste0("mkdir ", outdir_integ, "/report_images"))
  
  # create report
  rmarkdown::render(input = paste0(top_level_dir, "/bin/integration_module/integration_summary_report.Rmd"),
                    output_file = paste0(top_level_dir, "/", outdir_integ, "/summary_report.html"),
                    output_format = "html_document",
                    params = list(
                      directory = paste0(top_level_dir, "/", outdir_integ),
                      file = "peptide_info.tsv"
                    ))
  
  # zip all results files
  if (file.exists(paste0(outdir_integ, "/peptide_info.tsv")) && file.exists(paste0(outdir_integ, "/summary_report.html"))) {
    
    # create a zip file with results
    files_to_zip_int <- c("summary_report.html", "peptide_info.tsv", "report_images/",
                          "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                          "peptides.bed12", "ORFs.bed12", "transcripts.bed12","ncORF_stats.xlsx")
    
    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_int <- file.path("../integration_results.zip")
    
    # temp change the working dir to outdir_integ
    tmp_wd <- setwd(outdir_integ)
    
    # zip files
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
    
    # go back to starting dir
    setwd(top_level_dir)
    
  }
  
  
  
}

# main shiny app server
server <- function(input, output, session) {
  
  # store session ID
  # create session id tmp directory each time app is run
  session_id <- session$token
  print(paste0("Session: ", session_id))
  # create the dir
  system(paste0("mkdir ", session_id))
  
  # DATABASE MODULE
  
  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)
  
  # run database function when submit is pressed
  observeEvent(input$db_submit_button, {
    
    # ensure download button remains greyed out (if submit is re-pressed)
    shinyjs::disable("db_download_button")
    shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#d3d3d3';")
    # disable submit button after it is pressed
    session$sendCustomMessage("disableButton", list(id = "db_submit_button", spinnerId = "db-loading-container"))
    
    # run different servers depending on input type selected
    if (input$input_type == "fastq_input" & input$sequencing_type == "long-read") {
      fastq_server(input, output, session)
      bambu_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      bambu_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "fastq_input" & input$sequencing_type == "short-read") {
      fastq_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "bam_input" & input$sequencing_type == "short-read") {
      bam_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "gtf_input") {
      database_server(input, output, session)
    }
    
    # check if the zip file is created
    if (file.exists(paste0(session_id, "/database_results.zip"))) {
      file_available_db(TRUE)
    }
  })
  
  # enable download once files are available
  observe({
    if (file_available_db()) {
      shinyjs::enable("db_download_button")
      shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "db_submit_button", spinnerId = "db-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$db_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_database_results.zip")
    },
    content = function(file) {
      file.copy(paste0(session_id, "/database_results.zip"), file)
    }
  )
  
  # END DATABASE MODULE
  
  
  # PROTEOMICS MODULE
  
  # # create reactive value for the database zip
  # file_available_mm <- reactiveVal(FALSE)
  # 
  # 
  # # run database function when submit is pressed
  # observeEvent(input$proteomics_submit_button, { 
  #   session$sendCustomMessage("disableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # disable submit button
  #   proteomics_server(input, output, session)
  #   
  #   # check if the zip file is created
  #   if (file.exists("proteomics_output/proteomics_results.zip")) {
  #     file_available_mm(TRUE)
  #   }
  # })
  # 
  # # enable download once files are available
  # observe({
  #   if (file_available_mm()) {
  #     shinyjs::enable("proteomics_download_button")
  #     shinyjs::runjs("document.getElementById('proteomics_download_button').style.backgroundColor = '#4CAF50';")
  #     session$sendCustomMessage("enableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # re-enable submit button
  #   }
  # })
  # 
  # # download handler for the database results.zip file
  # output$proteomics_download_button <- downloadHandler(
  #   filename = function() {
  #     paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_proteomics_results.zip")
  #   },
  #   content = function(file) {
  #     file.copy("proteomics_output/proteomics_results.zip", file)
  #   }
  # )
  
  # END PROTEOMICS MODULE
  
  
  # INTEGRATION MODULE
  
  # create reactive value for the database zip
  file_available_integ <- reactiveVal(FALSE)
  
  # run integration function when submit is pressed
  observeEvent(input$integ_submit_button, { 
    
    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button
    
    # run integration server
    integration_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists(paste0(session_id, "/integration_results.zip"))) {
      file_available_integ(TRUE)
    }
    
  })
  
  # enable download once files are available
  observe({
    if (file_available_integ()) {
      shinyjs::enable("integ_download_button")
      shinyjs::runjs("document.getElementById('integ_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # re-enable submit button
    }
  })
  
  # download handler 
  output$integ_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_integration_results.zip")
    },
    content = function(file) {
      file.copy(paste0(session_id, "/integration_results.zip"), file)
    }
  )
  
  # END INTEGRATION MODULE
  
  
  # VISUALISATION MODULE
  
  data_storage <- reactiveValues()
  
  # function to import and process combined gtf
  import_and_preprocess_gtf <- function(gtf_path) {
    # import gtf and remove gene_id version number
    gtf_import <- rtracklayer::import(gtf_path, format = "gtf") %>%
      as_tibble() %>%
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    # separate gtf into multiple features
    list(
      gtf_import = gtf_import,
      res_tx_import = gtf_import %>% dplyr::filter(group_id == "transcripts"),
      res_ORF_import = gtf_import %>% dplyr::filter(group_id == "ORFs"),
      res_pep_import = gtf_import %>% dplyr::filter(group_id == "peptides")
    )
  }
  
  # function to update gene list
  update_gene_list <- function(res_tx_import, res_pep_import, uniq_map_peptides = FALSE, lncRNA_peptides = FALSE, novel_txs = FALSE, novel_txs_distinguished = FALSE, unann_orfs = FALSE, uorf_5 = FALSE, dorf_3 = FALSE) {
    
    # if check box is ticked
    if (uniq_map_peptides) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (lncRNA_peptides) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "lncRNA" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (novel_txs) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "novel" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (novel_txs_distinguished) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "novel" & peptide_ids_orf == TRUE & peptide_ids_transcript == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (unann_orfs) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(orf_type == "unannotated" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (uorf_5) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(localisation == "5UTR" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (dorf_3) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(localisation == "3UTR" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    # ensure gene_id overlaps
    ensembl_ids <- intersect(res_pep_import$gene_id, res_tx_import$gene_id)
    # filter for list
    genes_list <- res_tx_import %>% dplyr::filter(gene_id %in% ensembl_ids)
    
    # use gene name instead of ensembl
    if ("gene_name" %in% colnames(genes_list)) {
      unique(genes_list$gene_name)
    } else {
      unique(genes_list$gene_id)
    }
  }
  
  # function to update gene selection
  update_gene_selector <- function(session, genes_available) {
    updateSelectInput(session, "gene_selector", choices = genes_available)
  }
  
  # observe integration module file availability
  observe({
    # if the integration module has been run, use the output from this to display visualisation automatically
    if (file_available_integ()) {
      gtf_data <- import_and_preprocess_gtf(paste0(session_id, "/integ_output/combined_annotations.gtf"))
      data_storage$gtf_import <- gtf_data$gtf_import
      data_storage$res_tx_import <- gtf_data$res_tx_import
      data_storage$res_ORF_import <- gtf_data$res_ORF_import
      data_storage$res_pep_import <- gtf_data$res_pep_import
      
      genes_available <- update_gene_list(data_storage$res_tx_import, data_storage$res_pep_import, input$uniq_map_peptides, input$lncRNA_peptides, input$novel_txs, input$novel_txs_distinguished, input$unann_orfs, input$uorf_5, input$dorf_3)
      update_gene_selector(session, genes_available)
    }
  })
  
  # observe visualisation submit button
  observeEvent(input$vis_submit_button, {
    # if the submit button is pressed, import files
    session$sendCustomMessage("disableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container"))
    
    req(input$user_vis_gtf_file)
    
    gtf_data <- import_and_preprocess_gtf(input$user_vis_gtf_file$datapath)
    data_storage$gtf_import <- gtf_data$gtf_import
    data_storage$res_tx_import <- gtf_data$res_tx_import
    data_storage$res_ORF_import <- gtf_data$res_ORF_import
    data_storage$res_pep_import <- gtf_data$res_pep_import
    
    # if there are transcript counts and peptide intensities uploaded
    if (!is.null(input$user_vis_tx_count_file) & !is.null(input$user_pep_count_file)) {
      print("Counts detected")
      
      # import counts files
      data_storage$countst <- fread(input$user_vis_tx_count_file$datapath)
      data_storage$countsp <- fread(input$user_pep_count_file$datapath)
      
      # rename as per bambu counts output
      if ("TXNAME" %in% colnames(data_storage$countst) & "GENEID" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
        data_storage$countst$GENEID <- NULL
      } else if ("TXNAME" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
      }
      
      # filter GTF transcripts for those with counts
      data_storage$res_tx_import <- data_storage$res_tx_import %>% 
        dplyr::filter(transcript_id %in% data_storage$countst$transcript_id)
      
      # filter GTF ORFs for those with counts
      data_storage$res_ORF_import <- data_storage$res_ORF_import %>% 
        dplyr::filter(transcript_id %in% data_storage$countst$transcript_id)
      
      # match samples in both counts files
      sample_names <- intersect(colnames(data_storage$countsp), colnames(data_storage$countst))
      
      print("Samples with peptide intensities and transcript counts:")
      print(sample_names)
      
      # order sample names
      sample_names <- sample_names[order(match(sample_names, colnames(data_storage$countsp)))]
      
      # rename as per bambu counts output
      if ("Stripped.Sequence" %in% colnames(data_storage$countsp)) {
        data_storage$countsp$Peptide <- data_storage$countsp$Stripped.Sequence
        data_storage$countsp$Stripped.Sequence <- NULL
      } else if ("peptide" %in% colnames(data_storage$countsp)) {
        data_storage$countsp$Peptide <- data_storage$countsp$peptide
        data_storage$countsp$peptide <- NULL
      }
      
      # rarely, a peptide is in the data twice, so take max count value and get unique peptide IDs
      data_storage$countsp <- data_storage$countsp %>% 
        dplyr::select(Peptide, sample_names) %>% 
        dplyr::mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% 
        dplyr::group_by(Peptide) %>% 
        slice_max(sum) %>% dplyr::ungroup() %>% dplyr::select(-sum)
      
      # convert to matrix
      countsp_matrix <- as.matrix(data_storage$countsp[,-1])
      
      # set rownames as peptide IDs
      rownames(countsp_matrix) <- data_storage$countsp$Peptide
      
      # run VSN for normalisation of peptide intensities
      vsnp <- as.data.frame(justvsn(countsp_matrix))
      
      # create column of peptide IDs
      vsnp$peptide <- rownames(vsnp)
      
      # melt peptide intensities for plotting
      data_storage$countspm <- reshape2::melt(vsnp, id.vars = c("peptide"),
                                              variable.name = "sample_id", value.name = "count")
      
      # set levels for plotting
      data_storage$countspm$sample_id <- factor(as.character(data_storage$countspm$sample_id), level =  sample_names)
      
      # melt transcript counts for plotting
      data_storage$countstm <- reshape2::melt(data_storage$countst, id.vars = c("transcript_id"),
                                              variable.name = "sample_id", value.name = "count")
      # set levels for plotting
      data_storage$countstm$sample_id <- factor(as.character(data_storage$countstm$sample_id), level =  sample_names)
    }
    
    genes_available <- update_gene_list(data_storage$res_tx_import, data_storage$res_pep_import, input$uniq_map_peptides, input$lncRNA_peptides, input$novel_txs, input$unann_orfs, input$uorf_5, input$dorf_3)
    update_gene_selector(session, genes_available)
    
    # re-enable submit button after data is processed
    session$sendCustomMessage("enableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container"))
  })
  
  # make reactive list of checkbox filters
  checkbox_filters <- reactive({
    list(
      uniq_map_peptides = input$uniq_map_peptides,
      lncRNA_peptides = input$lncRNA_peptides,
      novel_txs = input$novel_txs,
      novel_txs_distinguished = input$novel_txs_distinguished,
      unann_orfs = input$unann_orfs,
      uorf_5 = input$uorf_5,
      dorf_3 = input$dorf_3
    )
  })
  
  # observe changes in any checkbox
  observe({
    # ensure GTFs are available
    req(data_storage$res_tx_import, data_storage$res_pep_import)
    
    # update gene list based on current checkbox states
    genes_available <- update_gene_list(
      data_storage$res_tx_import, 
      data_storage$res_pep_import, 
      uniq_map_peptides = checkbox_filters()$uniq_map_peptides,
      lncRNA_peptides = checkbox_filters()$lncRNA_peptides,
      novel_txs = checkbox_filters()$novel_txs,
      novel_txs_distinguished = checkbox_filters()$novel_txs_distinguished,
      unann_orfs = checkbox_filters()$unann_orfs,
      uorf_5 = checkbox_filters()$uorf_5,
      dorf_3 = checkbox_filters()$dorf_3
    )
    
    # update the gene selector with the new list
    update_gene_selector(session, genes_available)
  })
  
  # gene selection in drop down menu
  observeEvent(input$gene_selector, {
    
    req(input$gene_selector)
    
    # set gene name
    data_storage$gene_to_plot <- input$gene_selector
    print(data_storage$gene_to_plot)
    
    if (!is.null(input$user_vis_tx_count_file) & !is.null(input$user_pep_count_file)) { # if counts files are provided 
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import, data_storage$countstm, data_storage$countspm, min_intron_len=1000)
    } else { # if no counts
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import)
    }
    
    # print the plot
    output$plot <- renderPlot({
      suppressWarnings(print(data_storage$plot_obj))
    })
    
    # enable download button
    shinyjs::enable("vis_download_button")
    
  })
  
  # download handler for the plot
  output$vis_download_button <- downloadHandler(
    filename = function() {
      paste0(data_storage$gene_to_plot, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, h=10, w=20)
      print(data_storage$plot_obj)
      dev.off()
    }
  )
  
  # END VISUALISATION MODULE
  
  # remove session id tmp directory created each time app is run
  session$onSessionEnded(function() {
    if (dir.exists(session_id)) {
      unlink(session_id, recursive = TRUE)
    }
  })
  
}
