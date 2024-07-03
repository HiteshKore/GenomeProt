library(shiny)
library(shinyjs)

server <- function(input, output, session) {
  
  file_available_db <- reactiveVal(FALSE)
  file_available_integ <- reactiveVal(FALSE)
  data_storage <- reactiveValues()
  
  # database generation module
  observeEvent(input$db_submit_button, { 
    
    req(input$user_gtf_file)  # GTF is required
    
    gtf_path <- input$user_gtf_file$datapath
    
    if (!is.null(input$user_tx_count_file)) {
      tx_count_path <- input$user_tx_count_file$datapath
    } else {
      tx_count_path <- NULL
    }
    
    # run functions
    if (!is.null(tx_count_path)) {
      print("Counts")
      print(input$minimum_tx_count)
      filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, tx_counts=tx_count_path, min_count=input$minimum_tx_count)
    } else {
      print("No counts")
      filtered_gtf <- filter_custom_gtf(customgtf=gtf_path)
    }
    
    get_transcript_seqs("ORFome_transcripts.gtf", input$organism)
    # python cdhit
    
    # check files exist
    if (file.exists("ORFome_transcripts_nt.fasta") && file.exists("ORFome_transcripts.gtf")) {
      # create a zip file with results
      files_to_zip <- c("ORFome_transcripts_nt.fasta", "ORFome_transcripts.gtf")
      zipfile_path <- "results.zip"
      zip(zipfile = zipfile_path, files = files_to_zip)
      
      # check if the zip file is created
      if (file.exists(zipfile_path)) {
        file_available_db(TRUE)
      }
    }
  })
  
  observe({
    if (file_available_db()) {
      shinyjs::enable("db_download_button")
      shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#4CAF50';")
    }
  })
  
  # download handler for the results.zip file
  output$db_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_results.zip")
    },
    content = function(file) {
      file.copy("results.zip", file)
    }
  )
  #outputOptions(output, "dwonload_database_results ", suspendWhenHidden = FALSE)
  # end
  
  
  
  # integration module
  observeEvent(input$integ_submit_button, { 
    
    req(input$user_proteomics_file, input$user_post_gtf_file, input$user_fasta_file)  # GTF is required
    
    # file handling for different proteomics outputs
    #
    
    # Uncomment once fixed
    #system(paste0("Rscript map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -f ", input$user_fasta_file$datapath, " -g ", input$user_post_gtf_file$datapath))
    
    # check files exist
    if (file.exists("peptide_info.csv")) {
      # create a zip file with results
      files_to_zip_int <- c("peptide_info.csv", "peptides.gtf", "ORFs.gtf", "transcripts.gtf")
      zipfile_path_int <- "integration_results.zip"
      zip(zipfile = zipfile_path_int, files = files_to_zip_int)
      
      # check if the zip file is created
      if (file.exists(zipfile_path_int)) {
        file_available_integ(TRUE)
      }
    }
  })
  
  observe({
    if (file_available_integ()) {
      shinyjs::enable("integ_download_button")
      shinyjs::runjs("document.getElementById('integ_download_button').style.backgroundColor = '#4CAF50';")
    }
  })
  
  # download handler for the results.zip file
  output$downloadResultsIntegration <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_integration_results.zip")
    },
    content = function(file) {
      file.copy("integration_results.zip", file)
    }
  )
  
  # end
  
  
  
  # visualisation module
  observeEvent(input$vis_submit_button, { 
    
    req(input$user_tx_gtf_file, input$user_orf_gtf_file, input$user_pep_gtf_file)
    
    data_storage$res_tx_import <- rtracklayer::import(input$user_tx_gtf_file$datapath, format="gtf") %>% as_tibble() %>% 
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    
    data_storage$res_ORF_import <- rtracklayer::import(input$user_orf_gtf_file$datapath, format="gtf") %>% as_tibble()
    
    data_storage$res_pep_import <- rtracklayer::import(input$user_pep_gtf_file$datapath, format="gtf") %>% as_tibble() %>% 
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    
    if (!is.null(input$user_tx_count_file)) {
      
      print("Counts detected")
      data_storage$countst <- fread(input$user_tx_count_file$datapath)
      data_storage$countsp <- fread(input$user_pep_count_file$datapath)
      
      # when samples don't match
      sample_names <- intersect(colnames(data_storage$countsp), colnames(data_storage$countst))
      
      print("Samples with peptide intensities and transcript counts:")
      print(sample_names)
      
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
      
      data_storage$countstm <- reshape2::melt(data_storage$countst, id.vars = c("transcript_id"),
                                              variable.name = "sample_id", value.name = "count")
    } 
    
    # update genes available
    genes_available <- intersect(data_storage$res_pep_import$gene_id, data_storage$res_tx_import$gene_id)
    updateSelectInput(session, "gene_selector", choices = genes_available)
  })
  
  observeEvent(input$gene_selector, {
    
    req(input$gene_selector)
    
    data_storage$gene_to_plot <- input$gene_selector
    print(data_storage$gene_to_plot)
    
    if (!is.null(input$user_tx_count_file)) {
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
}
