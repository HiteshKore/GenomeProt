library(shiny)
library(shinyjs)

# internal server functions
fastq_server <- function(input, output, session) {
  
  req(input$user_reference_genome$datapath, input$user_reference_gtf$datapath, input$user_fastq_files$datapath)  # required
  
  outdir_bam="mapping_output"
  
  # check if the output directory exists
  if (dir.exists(outdir_bam)) {
    system(paste0("rm -rf ", outdir_bam, "/*"))
  } else {
    system(paste0("mkdir ", outdir_bam))
  }
  
  # create dataframe with sample details and add prefix column 
  user_fastq_files_df <- input$user_fastq_files %>% 
    mutate(file_prefix = str_replace_all(name, c("\\.fastq\\.gz$" = "", "\\.fastq$" = "", "\\.fq$" = "", "\\.fa$" = "", "\\.fasta$" = "")))
  
  print(user_fastq_files_df)
  
  # map reads
  if (input$sequencing_type == "long-read") {
  
    # generate index
    index_file <- paste0(outdir_bam, "/", str_replace_all(input$user_reference_genome$name, c("\\.fa$" = "", "\\.fasta$" = "")), ".mmi")
    
    minimap2_index_command <- paste0("minimap2 -ax splice:hq -d ", index_file, " ", input$user_reference_genome$datapath)
    #minimap2_index_command <- paste0("source activate py39; minimap2 -ax splice:hq -d ", index_file, " ", input$user_reference_genome$datapath)
    print(minimap2_index_command)
    system(minimap2_index_command)
  
    for (i in 1:nrow(user_fastq_files_df)) {
      fastq_file <- user_fastq_files_df$datapath[i]
      file_prefix <- user_fastq_files_df$file_prefix[i]
      
      minimap2_command <- paste0("minimap2 -t ", input$user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", index_file, " ", fastq_file, " | samtools view -bh -F 2308 | samtools sort -@ ", input$user_threads, " -o ", outdir_bam, "/", file_prefix, ".bam")
      #minimap2_command <- paste0("source activate py39; minimap2 -t ", input$user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", index_file, " ", fastq_file, " | samtools view -bh -F 2308 | samtools sort -@ ", input$user_threads, " -o ", outdir_bam, "/", file_prefix, ".bam")
      
      print(minimap2_command)
      system(minimap2_command)
    }

  } else if (input$sequencing_type == "short-read") {
    # short-read methods
  }
  
}

bambu_server <- function(input, output, session) {
  
  outdir_bam="mapping_output"
  
  if (input$input_type == "bam_input") {
    
    req(input$user_bam_files$datapath, input$user_reference_gtf$datapath, input$organism)  # required
    
    # create list of BAMs
    bam_file_list <- Rsamtools::BamFileList(as.vector(input$user_bam_files$datapath))
    print(bam_file_list)
    # get original names
    bam_file_names <- as.vector(input$user_bam_files$name)
    # remove bam extension
    bam_file_names <- str_remove(bam_file_names,".bam")
    # rename list to original names
    names(bam_file_list) <- bam_file_names
    print(bam_file_list)
    
  } else if (input$input_type == "fastq_input") {

    # create list of BAMs
    bam_files <- list.files(path = outdir_bam, "\\.bam$", full.names = TRUE)
    print(bam_files)
    
    bam_file_list <- Rsamtools::BamFileList(bam_files)
    print(bam_file_list)
    
    # remove bam extension
    bam_files_names <- list.files(path = outdir_bam, "\\.bam$", full.names = FALSE)
    bam_file_names <- str_remove(bam_files_names, ".bam")
    # rename list to original names
    names(bam_file_list) <- bam_file_names
    print(bam_file_list)
    
  }
  
  # run bambu function
  run_bambu_function(bam_file_list, input$user_reference_gtf$datapath, input$organism, input$user_threads)
  
  # run gffcompare
  
  #system(paste0("source activate py39; gffcompare -r ", input$user_reference_gtf$datapath, " bambu_output/bambu_transcript_annotations.gtf"))
  system(paste0("gffcompare -r ", input$user_reference_gtf$datapath, " bambu_output/bambu_transcript_annotations.gtf"))
  
  system(paste0("mv bambu_output/gffcmp.bambu_transcript_annotations.gtf.tmap bambu_output/gffcompare.tmap.txt"))
  system(paste0("rm gffcmp*"))
  
}

database_server <- function(input, output, session) {
  
  if (input$input_type == "gtf_input") {
    req(input$user_gtf_file, input$user_reference_gtf)  # GTFs required
    db_gtf_file <- input$user_gtf_file$datapath
    db_counts_file <- input$user_tx_count_file$datapath
  } else if ((input$input_type == "bam_input" | input$input_type == "fastq_input") & input$sequencing_type == "long-read") {
    db_gtf_file <- "bambu_output/bambu_transcript_annotations.gtf"
    db_counts_file <- "bambu_output/bambu_transcript_counts.txt"
  } else if ((input$input_type == "bam_input" | input$input_type == "fastq_input") & input$sequencing_type == "short-read") {
    # short-read methods
  }
  
  if (!is.null(db_counts_file)) {
    system(paste0("Rscript bin/database_module/generate_proteome.R -g ", db_gtf_file, " -r ", input$user_reference_gtf$datapath, " -c ", db_counts_file, 
                  " -m ", input$minimum_tx_count, " -o ", input$organism, " -l ", input$min_orf_length, " -u ", input$user_find_utr_5_orfs, " -d ", input$user_find_utr_3_orfs))
  } else {
    system(paste0("Rscript bin/database_module/generate_proteome.R -g ", db_gtf_file, " -r ", input$user_reference_gtf$datapath, 
                  " -o ", input$organism, " -l ", input$min_orf_length, " -u ", input$user_find_utr_5_orfs, " -d ", input$user_find_utr_3_orfs))
  }
  
  print("Generated ORFs")
  
  # set reference protein database
  if (input$organism == "human") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "mouse") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  }
  
  # run python script
  #system(paste0("source activate py39; python bin/database_module/annotate_proteome.py database_output/ref_transcripts_in_data.gtf ", ref_proteome, " database_output/ORFome_aa.txt database_output/proteome_database_transcripts.gtf database_output all"))
  system(paste0("python bin/database_module/annotate_proteome.py database_output/ref_transcripts_in_data.gtf ", ref_proteome, " database_output/ORFome_aa.txt database_output/proteome_database_transcripts.gtf database_output all"))
  
  print("Annotated proteome")
  
  # zip results
  if (file.exists("database_output/proteome_database.fasta") && file.exists("database_output/proteome_database_transcripts.gtf")) {
    if (input$input_type == "fastq_input") {
      bam_files <- list.files(path = "mapping_output/", "\\.bam$", full.names = TRUE)
      files_to_zip <- c(bam_files, "bambu_output/bambu_transcript_annotations.gtf", "bambu_output/bambu_transcript_counts.txt", "bambu_output/novel_transcript_classes.csv", "bambu_output/gffcompare.tmap.txt", "database_output/proteome_database.fasta", "database_output/proteome_database_metadata.txt", "database_output/proteome_database_transcripts.gtf")
    } else if (input$input_type == "bam_input") {
      files_to_zip <- c("bambu_output/bambu_transcript_annotations.gtf", "bambu_output/bambu_transcript_counts.txt", "bambu_output/novel_transcript_classes.csv", "bambu_output/gffcompare.tmap.txt", "database_output/proteome_database.fasta", "database_output/proteome_database_metadata.txt", "database_output/proteome_database_transcripts.gtf")
    } else if (input$input_type == "gtf_input") {
      files_to_zip <- c("database_output/proteome_database.fasta", "database_output/proteome_database_metadata.txt", "database_output/proteome_database_transcripts.gtf")
    }
    zipfile_path <- "database_output/database_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
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
  
  req(input$user_proteomics_file, input$user_post_gtf_file, input$user_fasta_file)  # GTF is required
  
  # run rscript
  system(paste0("Rscript bin/integration_module/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -f ", input$user_fasta_file$datapath, " -g ", input$user_post_gtf_file$datapath))
  
  # create report
  rmarkdown::render(input = "bin/integration_module/integration_summary_report.Rmd",
                    output_file = "../../integ_output/summary_report.html",
                    output_format = "html_document")
  
  # check files exist
  if (file.exists("integ_output/peptide_info.csv") && file.exists("integ_output/summary_report.html")) {
    # create a zip file with results
    files_to_zip_int <- c("integ_output/summary_report.html", "integ_output/peptide_info.csv", "integ_output/combined_annotations.gtf", "integ_output/peptides.bed12", "integ_output/ORFs.bed12", "integ_output/transcripts.bed12")
    zipfile_path_int <- "integ_output/integration_results.zip"
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
  }  
  
}

# shiny app server
server <- function(input, output, session) {
  
  # DATABASE MODULE
  
  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)
  
  # run database function when submit is pressed
  observeEvent(input$db_submit_button, {
    
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
      # short-read settings
    } else if (input$input_type == "bam_input" & input$sequencing_type == "short-read") {
      # short-read settings
    } else if (input$input_type == "gtf_input") {
      database_server(input, output, session)
    }
    
    # check if the zip file is created
    if (file.exists("database_output/database_results.zip")) {
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
      file.copy("database_output/database_results.zip", file)
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
  
  observeEvent(input$integ_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button
    integration_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists("integ_output/integration_results.zip")) {
      file_available_integ(TRUE)
    }
  })
  
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
      file.copy("integ_output/integration_results.zip", file)
    }
  )
  
  # END INTEGRATION MODULE
  
  
  # VISUALISATION MODULE
  data_storage <- reactiveValues()
  observeEvent(input$vis_submit_button, { 
    
    session$sendCustomMessage("disableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container")) # disable submit button
    
    req(input$user_vis_gtf_file)
    
    data_storage$gtf_import <- rtracklayer::import(input$user_vis_gtf_file$datapath, format="gtf") %>% as_tibble() %>% 
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    
    data_storage$res_tx_import <- data_storage$gtf_import %>% dplyr::filter(group_id == "transcripts")
    data_storage$res_ORF_import <- data_storage$gtf_import %>% dplyr::filter(group_id == "ORFs")
    data_storage$res_pep_import <- data_storage$gtf_import %>% dplyr::filter(group_id == "peptides")
    
    
    if (!is.null(input$user_vis_tx_count_file)) {
      
      print("Counts detected")
      # import counts
      data_storage$countst <- fread(input$user_vis_tx_count_file$datapath)
      data_storage$countsp <- fread(input$user_pep_count_file$datapath)
      
      # rename as per bambu counts output
      if ("TXNAME"  %in% colnames(data_storage$countst) & "GENEID" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
        data_storage$countst$GENEID <- NULL
      } else if ("TXNAME" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
      }
      
      # filter GTF transcripts for those with counts
      data_storage$res_tx_import <- data_storage$res_tx_import %>% 
        dplyr::filter(transcript_id %in% data_storage$countst$transcript_id)
      
      # match samples in both counts files
      sample_names <- intersect(colnames(data_storage$countsp), colnames(data_storage$countst))
      
      print("Samples with peptide intensities and transcript counts:")
      print(sample_names)
      
      sample_names <- sample_names[order(match(sample_names,colnames(data_storage$countsp)))]
      
      data_storage$countsp$Peptide <- data_storage$countsp$Stripped.Sequence
      
      # noted that sometimes a peptide is in the data twice, so take max count value
      data_storage$countsp <- data_storage$countsp %>% 
        dplyr::select(Peptide, sample_names) %>% 
        dplyr::mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% 
        dplyr::group_by(Peptide) %>% 
        slice_max(sum) %>% dplyr::ungroup() %>% dplyr::select(-sum)
      
      # VSN for normalisation 
      countsp_matrix <- as.matrix(data_storage$countsp[,-1])
      rownames(countsp_matrix) <- data_storage$countsp$Peptide
      
      # apply justvsn
      vsnp <- as.data.frame(justvsn(countsp_matrix))
      vsnp$peptide <- rownames(vsnp)
      
      # melt for plotting
      data_storage$countspm <- reshape2::melt(vsnp, id.vars = c("peptide"),
                                              variable.name = "sample_id", value.name = "count")
      
      # set levels for plotting
      data_storage$countspm$sample_id <- factor(as.character(data_storage$countspm$sample_id), level =  sample_names)
      
      data_storage$countstm <- reshape2::melt(data_storage$countst, id.vars = c("transcript_id"),
                                              variable.name = "sample_id", value.name = "count")
      # set levels for plotting
      data_storage$countstm$sample_id <- factor(as.character(data_storage$countstm$sample_id), level =  sample_names)
      
    } 
    
    # update genes available
    ensembl_ids <- intersect(data_storage$res_pep_import$gene_id, data_storage$res_tx_import$gene_id)
    
    genes_list <- data_storage$res_tx_import %>% 
      dplyr::filter(gene_id %in% ensembl_ids)
    
    genes_available <- unique(genes_list$gene_name)
    
    updateSelectInput(session, "gene_selector", choices = genes_available)
    
    session$sendCustomMessage("enableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container")) # re-enable submit button
    
  })
  
  observeEvent(input$gene_selector, {
    
    req(input$gene_selector)
    
    data_storage$gene_to_plot <- input$gene_selector
    print(data_storage$gene_to_plot)
    
    if (!is.null(input$user_vis_tx_count_file)) {
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import, data_storage$countstm, data_storage$countspm, min_intron_len=1000)
    } else {
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import)
    }
    
    # print the plot
    output$plot <- renderPlot({
      suppressWarnings(print(data_storage$plot_obj))
    })
    
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
  
}
