library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "GenomeProt",
                  dropdownMenu(type = "messages",
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/physiology/Parker-laboratory-Metabolic-Proteomics" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Parker Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Clark Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="mailto:genomeprot@outlook.com" target="_blank"><i class="fa fa-question"></i><h4>Support</h4><p>genomeprot@outlook.com</p></a></li>'))
                  )),
  # tabs
  dashboardSidebar(
    sidebarMenu(menuItem("Welcome", tabName = "welcome", icon = icon("house")),
                menuItem("Map FASTQs", tabName = "map_fastqs", icon = icon("dna")),
                menuItem("Identify isoforms", tabName = "run_bambu", icon = icon("magnifying-glass")),
                menuItem("Generate database", tabName = "db_generation", icon = icon("database")),
                menuItem("Analyse MS proteomics", tabName = "analyse_proteomics", icon = icon("gear")),
                menuItem("Integrate data", tabName = "integration", icon = icon("code-merge")),
                menuItem("Visualise results", tabName = "visualisation", icon = icon("chart-bar"))
    )
  ),
  # body
  dashboardBody(
    useShinyjs(),  # shinyjs
    tags$head(
      tags$style(HTML("
        .spinner {
          margin: 0 auto;
          width: 30px;
          height: 30px;
          border: 6px solid #ccc;
          border-top: 6px solid #333;
          border-radius: 50%;
          animation: spin 1s linear infinite;
        }
  
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
  
        .loading-container {
          display: none;
          text-align: center;
          margin-top: 20px;
        }
        
        #downloadResults {
          background-color: #4CAF50; /* Green */
          border: none;
          color: white;
          padding: 15px 32px;
          text-align: center;
          text-decoration: none;
          display: inline-block;
          font-size: 12px;
        }
        
        #downloadResults:disabled {
          background-color: #d3d3d3; /* Gray */
          color: #a9a9a9; /* Dark gray */
        }
        
        .spacing {
          margin-top: 20px;
        }
      ")),
      tags$script(HTML("
        Shiny.addCustomMessageHandler('disableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
          document.getElementById(params.spinnerId).style.display = 'block';
        });
  
        Shiny.addCustomMessageHandler('enableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
          document.getElementById(params.spinnerId).style.display = 'none';
        });
      "))
    ),
    tabItems(
      tabItem(tabName = "welcome",
              div(class = "jumbotron", style="background-image: url(dna-banner.svg); background-size: cover;", 
                  HTML("<center><h1>Welcome to GenomeProt</h1></center>"), 
                  HTML("<center><h3>An integrated proteogenomics data analysis platform</h3></center>"),
                  HTML("<center><h5>Developed by Hitesh Kore and Josie Gleeson at The University of Melbourne</h5></center>")
              ),
              fluidRow(
                column(12, 
                       div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:110%", 
                           div(class = "box-body", shiny::includeMarkdown("welcome-page-text.md"))
                       )
                )
              ),
              fluidRow(
                column(12, 
                       img(src = "images/LRPG-pipeline.png", width = "100%")
                )
              )
      ),
      tabItem(tabName = "map_fastqs", 
              h2("Map FASTQ files to the genome"),
              h5("NOTE: this step requires significant computation and time (>8 CPUs and high memory requirements)"),
              fluidRow(
                column(4,
                       selectInput("organism", label = "Organism:", 
                                   choices = list("human" = "human", "mouse" = "mouse"), 
                                   selected = "human"),
                       selectInput("sequencing_type", label = "Sequencing platform:", 
                                   choices = list("nanopore", "pacbio", "short-read"), 
                                   selected = "nanopore"),
                       fileInput("user_fastq_files", "Upload FASTQ file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE),
                       actionButton("map_fastqs_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("map_fastqs_download_button", "Download BAM file(s)", disabled = TRUE, style = "width:70%;")  # initially disabled
                )
              )
      ),
      tabItem(tabName = "run_bambu", 
              h2("Perform isoform discovery on BAM files with Bambu"),
              h5("NOTE: this step requires significant computation and time (>8 CPUs and high memory requirements)"),
              fluidRow(
                column(4,
                       #selectInput("organism", label = "Organism:", choices = list("human" = "human", "mouse" = "mouse"), selected = "human"),
                       fileInput("user_reference_genome", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_reference_gtf", "Upload reference GTF:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_bam_files", "Upload BAM file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE),
                       actionButton("bambu_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("bambu_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"),# initially disabled
                       div(id = "bambu-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "db_generation", 
              h2("Generate a custom proteogenomics database"),
              h5("Creates an amino acid FASTA of all ORFs in your data to use as input for MaxQuant/FragPipe etc."),
              fluidRow(
                column(4,
                       selectInput("organism", label = "Organism:", 
                                   choices = list("human" = "human", "mouse" = "mouse"), 
                                   selected = "human"),
                       #selectInput("startcodon", label = "Start codon:", choices = list("ATG" = "ATG", "ATG+CTG" = "ATG+CTG"), selected = "ATG"),
                       numericInput("min_orf_length", 
                                    label = "ORF Length (AA):", 
                                    value = 30),
                       checkboxInput("user_find_utr_orfs", label = "Search for 5' uORFs and 3' dORFs in UTRs of reference transcripts",
                                     value = FALSE, width = NULL),
                       fileInput("user_gtf_file", "Upload GTF file (from Bambu):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_ref_gtf_file", "Upload reference GTF file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_tx_count_file", "Upload transcript counts file (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       numericInput("minimum_tx_count", 
                                    label = "Minimum transcript expression threshold (leave blank if no counts file is uploaded, or for a default of 5):", 
                                    value = 5),
                       actionButton("db_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("db_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"),  # initially disabled
                       div(id = "db-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "analyse_proteomics", 
              h2("Run FragPipe with your custom proteogenomics database to analyse MS proteomics data"),
              h5("NOTE: this step requires significant computation and time (>8 CPUs and high memory requirements)"),
              fluidRow(
                column(4,
                       selectInput("organism", label = "Organism:", 
                                   choices = list("human" = "human", "mouse" = "mouse"), 
                                   selected = "human"),
                       numericInput("ORFLength", 
                                    label = "ORF Length (AA):", 
                                    value = 30),
                       fileInput("user_gtf_file", "Upload GTF file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_tx_count_file", "Upload transcript counts file (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       numericInput("minimum_tx_count", 
                                    label = "Minimum count per transcript (sum):", 
                                    value = 5),
                       actionButton("proteomics_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("proteomics_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "proteomics-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "integration", 
              h2("Integrate proteomics results with transcriptomics"),
              h5("Creates GTFs of peptides, ORFs and transcripts"),
              fluidRow(
                column(4,
                       fileInput("user_proteomics_file", "Upload proteomics file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_fasta_file", "Upload FASTA file used for proteomics:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_post_gtf_file", "Upload transcripts GTF used to build FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       actionButton("integ_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("integ_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "integ-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "visualisation", 
              h2("Visualise results"),
              h5("Plots your results using the GTFs created in the integration module."),
              fluidRow(
                column(4,
                       fileInput("user_tx_gtf_file", "Upload 'transcripts.gtf' file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_orf_gtf_file", "Upload 'ORFs.gtf' file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_pep_gtf_file", "Upload 'peptides.gtf' file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_vis_tx_count_file", "Upload transcript counts file (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_pep_count_file", "Upload peptide intensities file (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       actionButton("vis_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(8,
                       selectInput("gene_selector", "Select Gene", choices = NULL),
                       div(id = "vis-loading-container", class = "loading-container", div(class = "spinner")),
                       plotOutput("plot"),
                       downloadButton("vis_download_button", "Download plot", disabled = TRUE, class = "spacing")  # initially disabled
                )
              )
      )
    )
  ),  
  skin = "purple"
)
