library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  title = "GenomeProt",
  dashboardHeader(title = tags$img(
                    src = "images/GenomeProt_logo_Thach.png",
                    height = "40px",
                    style = "margin-left:1px",
                    style = "margin-right:1px"
                  ),
                  dropdownMenu(type = "messages",
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/physiology/Parker-laboratory-Metabolic-Proteomics" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Parker Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Clark Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="mailto:genomeprot@outlook.com" target="_blank"><i class="fa fa-question"></i><h4>Support</h4><p>genomeprot@outlook.com</p></a></li>'))
                  )),
  # tabs
  dashboardSidebar(width=200,
    sidebarMenu(menuItem("Welcome", tabName = "welcome", icon = icon("house")),
                menuItem("Generate database", tabName = "db_generation", icon = icon("database")),
                menuItem("Run proteomics analysis", tabName = "analyse_proteomics", icon = icon("search")),
                menuItem("Integrate data", tabName = "integration", icon = icon("code-merge")),
                menuItem("visualise results", tabName = "visualisation", icon = icon("eye")),
                menuItem("Quick help", tabName = "help", icon = icon("circle-question"))
    )
  ),
  # body
  dashboardBody(
    useShinyjs(),  # shinyjs
    tags$head(
      tags$style(HTML("
      .main-header .logo {
        width: 200px;
      }
      .main-header .navbar {
        margin-left: 200px;
      }
    ")),
      tags$link(rel = "shortcut icon", href = "favicon.ico"),
      tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "favicon.ico"),
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "/favicon-32x32.png"),
      tags$link(rel = "icon", type = "image/png", sizes = "16x16", href = "/favicon-16x16.png"),
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
        Shiny.addCustomMessageHandler('disableButtonOnly', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
        });

        Shiny.addCustomMessageHandler('disableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
          document.getElementById(params.spinnerId).style.display = 'block';
        });

        Shiny.addCustomMessageHandler('enableButtonOnly', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
        });

        Shiny.addCustomMessageHandler('enableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
          document.getElementById(params.spinnerId).style.display = 'none';
        });

        Shiny.addCustomMessageHandler('showStatusMessage', function(params) {
          var message = params.message;
          var color = params.color;
          var container = params.container;
          document.getElementById(container).innerText = message;
          document.getElementById(container).style.color = color;
        });

        Shiny.addCustomMessageHandler('clearStatusMessage', function(params) {
          var container = params.container;
          document.getElementById(container).innerText = '';
          document.getElementById(container).style.color = '';
        });

        var is_first_resize_ignored = false;
        window.addEventListener('message', handleMessage, false);

        function handleMessage(event)
        {
          let event_data = event.data;

          if (typeof(event_data) === 'string')
          {
            // tell IsoVis what tool it's embedded in
            if (event_data === 'Give tool name')
            {
              let isovis_window = document.getElementById('isovis_window');
              if (isovis_window && isovis_window.contentWindow)
                isovis_window.contentWindow.postMessage('Tool: GenomeProt', '*');
            }
            // set the height of the IsoVis iframe
            else if (event_data.startsWith('To resize: '))
            {
              if (!is_first_resize_ignored)
              {
                is_first_resize_ignored = true;
                return;
              }
              let new_window_height = Number.parseInt(event_data.substring('To resize: '.length));
              if (!Number.isNaN(new_window_height) && (new_window_height > 0))
              {
                if (new_window_height > window.innerHeight - 4)
                  new_window_height = window.innerHeight;
                else
                  new_window_height += 4;

                let isovis_window = document.getElementById('isovis_window');
                if (isovis_window)
                    isovis_window.style.height = `${new_window_height}px`;
              }
            }
          }
        }
      "))
    ),
    tabItems(
      tabItem(tabName = "welcome",
              fluidRow(
                column(12,
                       div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:110%",
                           div(class = "box-body", shiny::includeMarkdown("welcome-page-text.md")),
                           img(src = "images/workflow.png", width = "100%"),
                       )
                )
              )
      ),
      tabItem(tabName = "db_generation",
              h2("Generate a custom proteogenomics database"),
              h5("Creates an amino acid FASTA of all ORFs in your data to use as input for FragPipe/MaxQuant etc."),
              fluidRow(
                column(6,
                       # radio buttons for various options
                       radioButtons("sequencing_type", h5(tags$b("Select sequencing type:")),choices = c("Long-read (ONT, PacBio)" = "long-read", "Short-read" = "short-read")),
                       radioButtons("input_type", h5(tags$b("Select input type:")),choices = c("BAMs" = "bam_input","GTF (and/or transcript counts)" = "gtf_input")),
                       checkboxInput("vcf_option", "Incoporate SNVs into protein sequences", value = FALSE),

                       # organism
                       selectInput("organism", label = "Organism:",
                                   choices = list("Roundworm (C. elegans)" = "CAEEL", "Fruit fly (D. melanogaster)" = "DROME", "Human (H. sapiens)" = "HUMAN", "Mouse (M. musculus)" = "MOUSE", "Rat (R. rattus)" = "RAT", "Zebrafish (D. rerio)" = "DANRE"),
                                   selected = "HUMAN"),
                       # ORF type
                       selectInput("database_type", label = "ORFs to be included in proteomedb:", choices = list("canonical", "all"), selected = "all"),
                       # ORF length cutoff
                       numericInput("min_orf_length", label = "ORF length (amino acids):", value = 30),

                       # short ORF options
                       h5(tags$b("Find short (10 to 'ORF length' amino acids) ORFs in UTRs of reference transcripts:")),
                       checkboxInput("user_find_utr_5_orfs", label = "Upstream 5' ORFs",
                                     value = FALSE, width = NULL),
                       checkboxInput("user_find_utr_3_orfs", label = "Downstream 3' ORFs",
                                     value = FALSE, width = NULL),
                       numericInput("minimum_tx_count",
                                    label = "Minimum expression threshold (sum per transcript):",
                                    value = 5),

                       # Optional reference GTF file upload to override default
                       radioButtons("data_source","Reference annotation GTF:",
                                    choices = c( "Use preloaded GTF file" = "default","Upload Reference GTF file" = "user"), selected = "default"),

                       # GTF file upload appears only if user chooses to upload their own GTF file
                       conditionalPanel(condition = "input.data_source == 'user'",fileInput("reference_gtf_file","Upload GENCODE GTF File",buttonLabel = "Browse...",multiple = FALSE)),

                       # BAM input
                       conditionalPanel(
                         condition = "input.input_type == 'bam_input'",
                         numericInput("user_threads", label = "CPUs (Max 10):", value = 4, min = 1, max = 10, step = 1),

                         # short read
                         conditionalPanel(condition="(input.sequencing_type == 'short-read') & (input.vcf_option == true)",
                                          fileInput("user_reference_genome_bam", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE)),
                         conditionalPanel(condition="(input.sequencing_type == 'short-read')",
                                          fileInput("user_reference_genome_bam", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                                          fileInput("user_bam_files","Upload BAM Files", buttonLabel = "Browse...",multiple = TRUE)),

                         # radio button to load subsetted fasta file by default.
                         conditionalPanel(condition = "(input.sequencing_type == 'long-read') & (input.vcf_option == true)",
                           radioButtons("user_reference_genome_bam","Reference genome:",choices = c("Use preloaded reference genome file" = "default","Upload reference genome  FASTA" = "user"),selected = "default"),

                           # upload user specified genome fasta
                           conditionalPanel(condition = "(input.sequencing_type == 'long-read') & (input.user_reference_genome_bam == 'user')",
                             fileInput("user_reference_genome_bam","Upload reference genome  FASTA",buttonLabel = "Browse...",multiple = FALSE)),

                           # VCF file input conditionally shown
                           radioButtons("vcf_data","VCF file:",choices = c("Use preloaded VCF file" = "default","Upload VCF file" = "user"),selected = "default"),

                           # VCF file upload appears only if user chooses to upload their own VCF file
                           conditionalPanel(condition = "input.vcf_data == 'user'",
                             fileInput("user_vcf_file","Upload VCF File",buttonLabel = "Browse...", multiple = FALSE)),
                         ),

                         # radio button to load subsetted bam file by default
                         conditionalPanel(condition = "(input.sequencing_type == 'long-read')",
                           radioButtons("bam_data","BAM Files:",choices = c("Use preloaded BAM file" = "default","Upload BAM files" = "user"),selected = "default"),

                           # BAM file upload appears only if user chooses to upload their own BAM file
                           conditionalPanel(condition = "input.bam_data == 'user'",
                             fileInput("user_bam_files","Upload BAM Files", buttonLabel = "Browse...",multiple = TRUE))
                         ),
                       ),

                       # Custom GTF input (user generated GTF)
                       conditionalPanel(
                         condition = "input.input_type == 'gtf_input' & input.sequencing_type == 'long-read'",
                         fileInput("user_gtf_file", "Upload 'bambu_transcript_annotations.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                         fileInput("user_tx_count_file", "Upload 'bambu_transcript_counts.txt' (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.input_type == 'gtf_input' & input.sequencing_type == 'short-read'",
                         fileInput("user_tx_count_file", "Upload transcript counts:", NULL, buttonLabel = "Browse...", multiple = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.input_type == 'gtf_input' & input.vcf_option == true",
                         fileInput("user_genome_gtf", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE)
                       ),

                       actionButton("db_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("db_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "db-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "analyse_proteomics",
              h2("Perform proteomics searches using your custom database"),
              h5("Outputs: A zip file including 'peptide.tsv' (peptide search results), and if peptide quantification was enabled, 'report.pr_matrix.tsv' (peptide quantification results)."),
              h5("Note: This module only provides basic FragPipe settings. Users wishing to configure FragPipe beyond the level offered here or to use a proteomics pipeline other than FragPipe can perform their proteomics searches externally."),
              fluidRow(
                column(6,
                       h3("Perform a peptide search with FragPipe:"),
                       fileInput("fragpipe_prot_db_fasta_file", label = "Upload a proteome database FASTA file (e.g. 'proteome_database.fasta'):", buttonLabel = "Browse...", multiple = FALSE, accept = c(".fasta")),
                       checkboxInput("user_add_contaminants", label = "Add contaminants into the proteome database?",
                                     value = TRUE, width = NULL),
                       checkboxInput("user_perform_quantification", label = "Perform peptide quantification after the peptide search?",
                                     value = FALSE, width = NULL),
                       h5(tags$b("Upload mass spectrometry data files and select their data types:")),
                       actionButton("add_mass_spec_file_button", "+ Add mass spectrometry data file", class = "btn btn-warning"),
                       actionButton("remove_mass_spec_file_button", "- Remove mass spectrometry data file", class = "btn"),
                       div(id = "mass_spec_file_list"),
                       selectInput("protease1", label = "Protease 1:",
                                   choices = list("stricttrypsin (cuts KR, sense C)" = "stricttrypsin",
                                                  "trypsin (cuts KR, no cuts P, sense C)" = "trypsin",
                                                  "trypsin_gluc (cuts DEKR, no cuts P, sense C)" = "trypsin_gluc",
                                                  "gluc (cuts DE, no cuts P, sense C)" = "gluc",
                                                  "lysc (cuts K, no cuts P, sense C)" = "lysc",
                                                  "lysn (cuts K, sense N)" = "lysn",
                                                  "argc (cuts R, no cuts P, sense C)" = "argc",
                                                  "aspn (cuts D, sense N)" = "aspn"),
                                   selected = "stricttrypsin (cuts KR, sense C)"),
                       selectInput("protease2", label = "Protease 2 (optional; must be different from protease 1):",
                                   choices = list("none" = "none",
                                                  "stricttrypsin (cuts KR, sense C)" = "stricttrypsin",
                                                  "trypsin (cuts KR, no cuts P, sense C)" = "trypsin",
                                                  "trypsin_gluc (cuts DEKR, no cuts P, sense C)" = "trypsin_gluc",
                                                  "gluc (cuts DE, no cuts P, sense C)" = "gluc",
                                                  "lysc (cuts K, no cuts P, sense C)" = "lysc",
                                                  "lysn (cuts K, sense N)" = "lysn",
                                                  "argc (cuts R, no cuts P, sense C)" = "argc",
                                                  "aspn (cuts D, sense N)" = "aspn"),
                                   selected = "none"),
                       numericInput("fragpipe_cpu_threads", label = "Number of CPU threads to use (specify '0' for FragPipe to use <number of cores in system - 1> threads)", value = 1, step = 1),
                       numericInput("fragpipe_memory_limit", label = "Memory limit in GB to use (specify '0' to let FragPipe decide)", value = 15, step = 1),
                       actionButton("fragpipe_submit_button", "Run FragPipe", class = "btn btn-info")
                ),
                column(6,
                       h3("Download FragPipe results:"),
                       downloadButton("fragpipe_download_button", "Download FragPipe results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "fragpipe-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "fragpipe-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "integration",
              h2("Integrate proteomics results with transcriptomics"),
              h5("Creates BED12s and GTFs of peptides, ORFs and transcripts for visualisation and produces summary data"),
              fluidRow(
                column(6,
                       fileInput("user_proteomics_file", "Upload proteomics results:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_metadata_file", "Upload 'proteome_database_metadata.txt':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_post_gtf_file", "Upload 'proteome_database_transcripts.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE),
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
              h2("Visualise results with IsoVis"),
              h5("The IsoVis website is displayed below for convenience. It is also accessible directly at: https://isomix.org/isovis/"),
              h5(actionLink("show_isovis_steps", "Instructions for using IsoVis")),
              conditionalPanel(
                condition = "input.show_isovis_steps % 2 == 1",
                p("Step 1: Click 'Upload data'."),
                p("Step 2: For the 'Transcript data' file, upload 'combined_annotations.gtf'. For the 'Transcript counts' file, upload 'bambu_transcript_counts.txt'."),
                p("Step 3: For the 'Peptide intensities' file, upload the peptide intensities file from the proteomics pipeline you used (e.g. 'report.pr_matrix.tsv'), then click 'Apply'."),
                p("Step 4: Type the symbol or ID of a gene to view, select it from the list of results displayed, then either press enter or click '>'."),
                p("Step 5: To see the mappings of peptides to open reading frames, click on the 'Stack options' dropdown menu and select 'Peptide mapping'.")
              ),
              fluidRow(
                column(12,
                    tags$iframe(id = "isovis_window",
                                src = "https://isomix.org/isovis/",
                                width = "100%",
                                height = "950px",
                                style = "border:none;"))
              )
      ),
      tabItem(
        tabName = "help",
        fluidRow(
          column(
            width = 12,
            box(
              width = 20,
              status = "primary",
              solidHeader = TRUE,
              tags$iframe(src = "GenomeProt_help.html",
                          width = "100%",
                          height = "950px",
                          style = "border:none;")
            )
          )
        )
      )
    )
  ),
  skin = "purple"
)