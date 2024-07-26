library(shiny)
library(shinyjs)

# internal server functions
fastq_server <- function(input, output, session) {
  
  req(input$user_threads,input$user_reference_genome$datapath,input$user_fastq_files$datapath)  # required
  
  outdir_bam="fastq_output/"
  # check if the output directory exists
  if (dir.exists(outdir_bam)) {
    system(paste0("rm -rf ", outdir_bam, "*"))
  } else {
    system(paste0("mkdir ", outdir_bam))
  }
  
  # create dataframe with sample details and add prefix column 
  user_fastq_files_df <- input$user_fastq_files %>% 
    mutate(file_prefix = str_replace_all(name, c("\\.fastq\\.gz$" = "", "\\.fastq$" = "")))
  
  # map reads
  if (input$sequencing_type=="long-read") {
  
    # generate index
    index_file <- paste0("fastq_output/", str_replace_all(input$user_reference_genome$name, c("\\.fa$" = "", "\\.fasta$" = "")), ".mmi")
    
    #minimap2_index_command <- paste0("minimap2 -ax splice:hq -d ", index_file, " ", input$user_reference_genome$datapath)
    minimap2_index_command <- paste0("source activate IsoLamp; minimap2 -ax splice:hq -d ", index_file, " ", input$user_reference_genome$datapath)
    print(minimap2_index_command)
    system(minimap2_index_command)
  
    for (i in 1:nrow(user_fastq_files_df)) {
      fastq_file <- user_fastq_files_df$datapath[i]
      file_prefix <- user_fastq_files_df$file_prefix[i]
      
      #minimap2_command <- paste0("minimap2 -t ", input$user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", index_file, " ", fastq_file, " | samtools view -bh -F 2308 | samtools sort -@ ", input$user_threads, " -o ", outdir_bam, file_prefix, ".bam")
      minimap2_command <- paste0("source activate IsoLamp; minimap2 -t ", input$user_threads, " -ax splice:hq --sam-hit-only --secondary=no ", index_file, " ", fastq_file, " | samtools view -bh -F 2308 | samtools sort -@ ", input$user_threads, " -o ", outdir_bam, file_prefix, ".bam")
      
      print(minimap2_command)
      system(minimap2_command)
    }

  }
  
  # check bam files exist
  bam_files <- list.files(path = outdir_bam, "\\.bam$", full.names = TRUE)
  
  if (length(bam_files)>0) {
    # create a zip file with results
    zipfile_path <- paste0(outdir_bam, "bam_results.zip")
    zip(zipfile = zipfile_path, files = bam_files)
  }
  
}

bambu_server <- function(input, output, session) {
  
  req(input$user_bam_files, input$user_reference_gtf$datapath, input$bambu_organism)  # required
  
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
  
  # run bambu function
  run_bambu_function(bam_file_list, input$user_reference_gtf$datapath, input$organism)
  
  # run gffcompare
  
  system(paste0("source activate IsoLamp; gffcompare -r ", input$user_reference_gtf$datapath, " bambu_output/bambu_transcript_annotations.gtf"))
  #system(paste0("gffcompare -r ", input$user_reference_genome$datapath, " bambu_output/bambu_transcript_annotations.gtf"))
  
  system(paste0("mv bambu_output/gffcmp.bambu_transcript_annotations.gtf.tmap bambu_output/gffcompare.tmap.txt"))
  system(paste0("rm gffcmp*"))
  
  # check files exist
  if (file.exists("bambu_output/bambu_transcript_annotations.gtf") && file.exists("bambu_output/gffcompare.tmap.txt")) {
    # create a zip file with results
    files_to_zip <- c("bambu_output/bambu_transcript_annotations.gtf", "bambu_output/bambu_transcript_counts.txt", "bambu_output/novel_transcript_classes.csv", "bambu_output/gffcompare.tmap.txt")
    zipfile_path <- "bambu_output/bambu_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
  }
  
  # if short-reads are supplied this won't work
  # short-read data should skip isoform identification and just perform quantification, then filter based on counts
  
}

database_server <- function(input, output, session) {
  
  req(input$user_gtf_file, input$user_ref_gtf_file)  # GTFs required
  
  if (!is.null(input$user_tx_count_file)) {
    system(paste0("Rscript bin/database_module/generate_proteome.R -g ", input$user_gtf_file$datapath, " -r ", input$user_ref_gtf_file$datapath, " -c ", input$user_tx_count_file$datapath, 
                  " -m ", input$minimum_tx_count, " -o ", input$organism, " -l ", input$min_orf_length, " -u ", input$user_find_utr_orfs))
  } else {
    system(paste0("Rscript bin/database_module/generate_proteome.R -g ", input$user_gtf_file$datapath, " -r ", input$user_ref_gtf_file$datapath, 
                  " -o ", input$organism, " -l ", input$min_orf_length, " -u ", input$user_find_utr_orfs))
  }
  
  print("Generated ORFs")
  
  # set reference protein database
  if (input$organism == "human") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "mouse") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  }
  
  # run python script
  system(paste0("source activate py39; python bin/database_module/annotate_proteome.py db_output/ref_transcripts_in_data.gtf ", ref_proteome, " db_output/ORFome_aa.txt db_output/proteome_database_transcripts.gtf db_output all"))
  #system(paste0("python bin/database_module/annotate_proteome.py db_output/ref_transcripts_in_data.gtf ", ref_proteome, " db_output/ORFome_aa.txt db_output/proteome_database_transcripts.gtf db_output all"))
  
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
  
  system(paste0("source activate mm_env; metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  #system(paste0("metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  
  # check files exist
  if (file.exists("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv") && file.exists("proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")) {
    # create a zip file with results
    files_to_zip <- c("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv", "proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
    zipfile_path <- "proteomics_output/proteomics_results.zip"
    zip(zipfile = zipfile_path, files = files_to_zip)
  }
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
  
  # MAP FASTQ MODULE
  
  # check if the zip file is created
  file_available_bam <- reactiveVal(FALSE)
  
  observeEvent(input$map_fastqs_submit_button, { 
    
    session$sendCustomMessage("disableButton", list(id = "map_fastqs_submit_button", spinnerId = "fastq-loading-container")) # disable submit button
    fastq_server(input, output, session)
    if (file.exists("fastq_output/bam_results.zip")) {
      file_available_bam(TRUE)
    }
  })
  
  # enable download once files are available
  observe({
    if (file_available_bam()) {
      shinyjs::enable("map_fastqs_download_button")
      shinyjs::runjs("document.getElementById('map_fastqs_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "map_fastqs_submit_button", spinnerId = "fastq-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$map_fastqs_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_bam_results.zip")
    },
    content = function(file) {
      file.copy("fastq_output/bam_results.zip", file)
    }
  )
  
  # END FASTQ MODULE
  
  
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
    genes_available <- intersect(data_storage$res_pep_import$gene_id, data_storage$res_tx_import$gene_id)
    
    # to use names
    # library(biomaRt)
    # 
    # mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
    # 
    # gene_names <- getBM(filters= "ensembl_gene_id",
    #                          attributes= c("ensembl_gene_id","hgnc_symbol"),
    #                          values = genes_available,
    #                          mart = mart)
    # 
    # data_storage$genes_available <- gene_names$ensembl_gene_id
    # names(data_storage$genes_available) <- gene_names$hgnc_symbol
    
    # change back from names?
    #updateSelectInput(session, "gene_selector", choices = names(data_storage$genes_available))
    
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
