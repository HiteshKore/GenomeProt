library(shiny)
library(shinyjs)

# internal server functions
fastq_server <- function(input, output, session) {
  
  
  renderTable(input$user_fastq_files)
  
  # call minimap2 script and wait for BAM output
  user_threads <- 12
  user_genome <- "genome.fa" 
  fastq_file <- "test"
  
  system(paste0("minimap2 -t ", user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", user_genome, " ", fastq_file, ".fastq | samtools view -bh -F 2308 | samtools sort -@ ", user_threads, " -o ", fastq_file, ".bam"))
  #system(paste0("Rscript bin/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -f ", input$user_fasta_file$datapath, " -g ", input$user_post_gtf_file$datapath))
  
  # short-reads, call STAR

}

bambu_server <- function(input, output, session) {
  
  req(input$user_bam_files, input$user_reference_gtf$datapath, input$user_reference_genome$datapath)  # required

  library(bambu)
  # create list of BAMs
  bam_file_list <- Rsamtools::BamFileList(as.vector(input$user_bam_files$datapath))
  # get original names
  bame_file_names <- as.vector(input$user_bam_files$name)
  # remove bam extension
  bame_file_names <- str_remove(bame_file_names,".bam")
  # rename list to original names
  names(bam_file_list) <- bame_file_names
  
  # run bambu function
  run_bambu_function(bam_file_list, input$user_reference_gtf$datapath, input$user_reference_genome$datapath)
  
  # check files exist
  if (file.exists("bambu_output/bambu_transcript_annotations.gtf")) {
    # create a zip file with results
    files_to_zip <- c("bambu_output/bambu_transcript_annotations.gtf", "bambu_output/bambu_transcript_counts.txt", "bambu_output/novel_transcript_classes.csv")
    zipfile_path <- "bambu_output/bambu_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
  }
  
  # if short-reads are supplied this won't work
  # short-read data should skip isoform identification and just perform quantification, then filter based on counts

}

database_server <- function(input, output, session) {
  
  req(input$user_gtf_file, input$user_ref_gtf_file)  # GTF is required
  
  gtf_path <- input$user_gtf_file$datapath
  reference_gtf <- input$user_ref_gtf_file$datapath
  
  system("mkdir db_output")
  
  if (!is.null(input$user_tx_count_file)) {
    tx_count_path <- input$user_tx_count_file$datapath
  } else {
    tx_count_path <- NULL
  }
  
  # run filter_custom_gtf, if counts are present, supply them
  if (!is.null(tx_count_path)) {
    filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, referencegtf=reference_gtf, tx_counts=tx_count_path, min_count=input$minimum_tx_count)
  } else {
    filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, referencegtf=reference_gtf)
  }
  
  # currently calls a function defined in 'functions.R'
  # run filter_custom_gtf, if counts are present, supply them
  get_transcript_seqs(filteredgtf="db_output/proteome_database_transcripts.gtf", organism=input$organism, orf_len=input$min_orf_length, find_UTR_orfs=input$user_find_utr_orfs, referencegtf=reference_gtf)
  
  if (input$organism == "human") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "mouse") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  }
  
  gtf_type <- "GENCODE"
  
  # run python script
  #system(paste0("source activate py39; python bin/annotate_proteome.py db_output/ref_transcripts_in_data.gtf ", ref_proteome, " db_output/ORFome_aa.txt db_output/proteome_database_transcripts.gtf ", gtf_type))
  system(paste0("python bin/annotate_proteome.py db_output/ref_transcripts_in_data.gtf ", ref_proteome, " db_output/ORFome_aa.txt db_output/proteome_database_transcripts.gtf ", gtf_type))
  print("Annotated proteome")
  
  # check files exist
  if (file.exists("db_output/proteome_database.fasta") && file.exists("db_output/proteome_database_transcripts.gtf")) {
    # create a zip file with results
    files_to_zip <- c("db_output/proteome_database.fasta", "db_output/proteome_database_metadata.txt", "db_output/proteome_database_transcripts.gtf")
    zipfile_path <- "db_output/database_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
  }
  
}

proteomics_server <- function(input, output, session) {
  
  req(input$user_mm_data, input$user_mm_fasta)
  
  #mass_spec_data <- as.vector(input$user_mm_data$datapath)
  #print(mass_spec_data)
  #print(input$user_mm_data$name)
  
  # get directory path
  dir_path <- dirname(input$user_mm_data$datapath[1])
  print(dir_path)
  
  mass_spec_names <- c(input$user_mm_data$name)
  print(mass_spec_names)
  
  # new file paths
  new_file_paths <- file.path(dir_path, mass_spec_names)

  # rename files
  file.rename(input$user_mm_data$datapath, new_file_paths)
  
  # copy config files
  # replace all lines:
  # 'MaxThreadsToUsePerFile = 3'
  # with user set threads
  
  #system(paste0("source activate mm_env; metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  system(paste0("metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  
  # check files exist
  if (file.exists("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv") && file.exists("proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")) {
    # create a zip file with results
    files_to_zip <- c("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv", "proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
    zipfile_path <- "proteomics_output/proteomics_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
  }
}

integration_server <- function(input, output, session) {
  # NOTE: create results dir
  req(input$user_proteomics_file, input$user_post_gtf_file, input$user_fasta_file)  # GTF is required
  system("mkdir integ_output")
  # file handling for different proteomics pipelines
  
  # un-comment once fixed
  system(paste0("Rscript bin/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -f ", input$user_fasta_file$datapath, " -g ", input$user_post_gtf_file$datapath))
  
  # check files exist
  if (file.exists("integ_output/peptide_info.csv")) {
    # create a zip file with results
    files_to_zip_int <- c("integ_output/peptide_info.csv", "integ_output/peptides.gtf", "integ_output/ORFs.gtf", "integ_output/transcripts.gtf")
    zipfile_path_int <- "integration_results.zip"
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
  }  
  
}

visualisation_server <- function(input, output, session) {
  
  # unsure if the following will break the code if run in sep function
  # data_storage <- reactiveValues()
  
}

# shiny app server
server <- function(input, output, session) {
  
  # IDENTIFY ISOFORMS MODULE
  file_available_bambu <- reactiveVal(FALSE)
  
  observeEvent(input$bambu_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "bambu_submit_button", spinnerId = "bambu-loading-container")) # disable submit button
    bambu_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists("bambu_output/bambu_results.zip")) {
      file_available_bambu(TRUE)
      
    }
  })
  
  # enable download once files are available
  observe({
    if (file_available_bambu()) {
      shinyjs::enable("bambu_download_button")
      shinyjs::runjs("document.getElementById('bambu_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "bambu_submit_button", spinnerId = "bambu-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$bambu_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_bambu_results.zip")
    },
    content = function(file) {
      file.copy("bambu_output/bambu_results.zip", file)
    }
  )
  
  file_available_bambu <- reactiveVal(FALSE)
  
  observeEvent(input$bambu_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "bambu_submit_button", spinnerId = "bambu-loading-container")) # disable submit button
    bambu_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists("bambu_output/bambu_results.zip")) {
      file_available_bambu(TRUE)
      
    }
  })
  
  # enable download once files are available
  observe({
    if (file_available_bambu()) {
      shinyjs::enable("bambu_download_button")
      shinyjs::runjs("document.getElementById('bambu_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "bambu_submit_button", spinnerId = "bambu-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$bambu_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_bambu_results.zip")
    },
    content = function(file) {
      file.copy("bambu_output/bambu_results.zip", file)
    }
  )
  # END IDENTIFY MODULE
  
  
  # PROTEOMICS MODULE
  
  # create reactive value for the database zip
  file_available_mm <- reactiveVal(FALSE)
  
  # run database function when submit is pressed
  observeEvent(input$proteomics_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # disable submit button
    proteomics_server(input, output, session)

    # check if the zip file is created
    if (file.exists("proteomics_output/proteomics_results.zip")) {
      file_available_mm(TRUE)
    }
  })
  
  # enable download once files are available
  observe({
    if (file_available_mm()) {
      shinyjs::enable("proteomics_download_button")
      shinyjs::runjs("document.getElementById('proteomics_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$proteomics_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_proteomics_results.zip")
    },
    content = function(file) {
      file.copy("proteomics_output/proteomics_results.zip", file)
    }
  )
  
  # END PROTEOMICS MODULE
  
  
  # DATABASE MODULE
  
  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)
  
  # run database function when submit is pressed
  observeEvent(input$db_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "db_submit_button", spinnerId = "db-loading-container")) # disable submit button
    database_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists("db_output/database_results.zip")) {
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
      file.copy("db_output/database_results.zip", file)
    }
  )
  
  # END DATABASE MODULE
  
  
  # INTEGRATION MODULE
  
  # create reactive value for the database zip
  file_available_integ <- reactiveVal(FALSE)
  
  observeEvent(input$integ_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button
    integration_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists("integration_results.zip")) {
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
  output$downloadResultsIntegration <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_integration_results.zip")
    },
    content = function(file) {
      file.copy("integration_results.zip", file)
    }
  )
  
  # END INTEGRATION MODULE
  
  
  # VISUALISATION MODULE
  data_storage <- reactiveValues()
  observeEvent(input$vis_submit_button, { 
    session$sendCustomMessage("disableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container")) # disable submit button
    req(input$user_tx_gtf_file, input$user_orf_gtf_file, input$user_pep_gtf_file)
    
    data_storage$res_tx_import <- rtracklayer::import(input$user_tx_gtf_file$datapath, format="gtf") %>% as_tibble() %>% 
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    
    data_storage$res_ORF_import <- rtracklayer::import(input$user_orf_gtf_file$datapath, format="gtf") %>% as_tibble()
    
    data_storage$res_pep_import <- rtracklayer::import(input$user_pep_gtf_file$datapath, format="gtf") %>% as_tibble() %>% 
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    
    if (!is.null(input$user_vis_tx_count_file)) {
      
      print("Counts detected")
      data_storage$countst <- fread(input$user_vis_tx_count_file$datapath)
      data_storage$countsp <- fread(input$user_pep_count_file$datapath)
      
      # when samples don't match
      sample_names <- intersect(colnames(data_storage$countsp), colnames(data_storage$countst))
      #sample_names <- intersect(colnames(countsp), colnames(countst))

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
      
      # VSN
      countsp_matrix <- as.matrix(data_storage$countsp[,-1])
      rownames(countsp_matrix) <- data_storage$countsp$Peptide
      
      # apply justvsn
      vsnp <- as.data.frame(justvsn(countsp_matrix))
      vsnp$peptide <- rownames(vsnp)
      
      # melt for plotting
      data_storage$countspm <- reshape2::melt(vsnp, id.vars = c("peptide"),
                                              variable.name = "sample_id", value.name = "count")
      
      data_storage$countspm$sample_id <- factor(as.character(data_storage$countspm$sample_id), level =  sample_names)
      
      data_storage$countstm <- reshape2::melt(data_storage$countst, id.vars = c("transcript_id"),
                                              variable.name = "sample_id", value.name = "count")
      
      data_storage$countstm$sample_id <- factor(as.character(data_storage$countstm$sample_id), level =  sample_names)
      
    } 
    
    # update genes available
    genes_available <- intersect(data_storage$res_pep_import$gene_id, data_storage$res_tx_import$gene_id)
    updateSelectInput(session, "gene_selector", choices = genes_available)
    session$sendCustomMessage("enableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container")) # re-enable submit button
  })
  
  observeEvent(input$gene_selector, {
    
    req(input$gene_selector)
    
    data_storage$gene_to_plot <- input$gene_selector
    print(data_storage$gene_to_plot)
    
    if (!is.null(input$user_vis_tx_count_file)) {
      print("counts")
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import, data_storage$countstm, data_storage$countspm, min_intron_len=1000)
    } else {
      print("no counts")
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
      pdf(file)
      print(data_storage$plot_obj)
      dev.off()
    }
  )
  # END VISUALISATION MODULE
  
}
