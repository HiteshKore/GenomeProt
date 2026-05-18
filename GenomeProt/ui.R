library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = tags$img( ##title = "GenomeProt"
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
                # commenting out proteomics for now
                #menuItem("Analyse MS proteomics", tabName = "analyse_proteomics", icon = icon("gear")),
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
                       radioButtons("input_type", h5(tags$b("Select input type:")),choices = c( #"FASTQs" = "fastq_input",
                         "BAMs" = "bam_input","GTF (and/or transcript counts)" = "gtf_input")),
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

                       # conditionl panel for FASTQ input
                       conditionalPanel(condition = "input.input_type == 'fastq_input'",
                         numericInput("user_threads", label = "CPUs (Max 10):", value = 4, min = 1, max = 10, step = 1),
                         #h5("Map FASTQs, identify (in long-reads) and quantify isoforms, and generate the database"),
                         fileInput("user_reference_genome", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                         conditionalPanel(condition = "input.sequencing_type == 'short-read'",
                           fileInput("transcriptome_file", "Upload reference transcriptome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE)
                         ),
                         fileInput("user_fastq_files", "Upload FASTQ file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE)
                       ),

                       # BAM input
                       conditionalPanel(
                         condition = "input.input_type == 'bam_input'",
                         numericInput("user_threads", label = "CPUs (Max 10):", value = 4, min = 1, max = 46, step = 1),

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
      # tabItem(tabName = "analyse_proteomics",
      #         h2("Run MetaMorpheus with your custom proteogenomics database to analyse MS proteomics data"),
      #         h5("NOTE: this step requires significant computation and time (>8 CPUs and high memory requirements)"),
      #         fluidRow(
      #           column(4,
      #                  selectInput("protease", label = "Protease:",
      #                              choices = list("trypsin" = "trypsin"),
      #                              selected = "trypsin"),
      #                  numericInput("mm_cpu",
      #                               label = "CPUs",
      #                               value = 1),
      #                  fileInput("user_mm_fasta", "Upload 'proteome_database.fasta'", NULL, buttonLabel = "Browse...", multiple = FALSE),
      #                  fileInput("user_mm_data", "Upload mzML/raw file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE),
      #                  actionButton("proteomics_submit_button", "Submit", class = "btn btn-primary")
      #           ),
      #           column(6,
      #                  HTML("<h3>Download your results:</h3>"),
      #                  downloadButton("proteomics_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
      #                  div(id = "proteomics-loading-container", class = "loading-container", div(class = "spinner"))
      #           )
      #         )
      # ),
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
                p("Step 1: Click 'Upload data'. For the 'Stack data' upload 'transcripts_and_ORFs_for_isovis.gtf'. For the 'Heatmap data' upload 'bambu_transcript_counts.txt'."),
                p("Step 2: Check the box 'Show peptide data upload options'."),
                p("Step 3: For the 'Peptide sites data' upload 'peptides.bed12'. For the 'Peptide intensities data' upload the intensities file. Then click 'Apply'."),
                p("Step 4: Type a gene to view and press enter."),
                p("Step 5: Click 'Stack options' and select 'Peptide sites' from the drop-down menu.")
              ),
              fluidRow(
                column(12,
                    tags$iframe(src = "https://isomix.org/isovis/",
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
              includeHTML("GenomeProt_help.html")
              #tags$iframe(src="GenomeProt_help.html",
              #width= "100%",
              #height="750px",
              #style="border:none;"
              #)
            )
          )
        )
      )
    )
  ),
  skin = "purple"
)