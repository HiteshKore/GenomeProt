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
  dashboardSidebar(width = 200,
    sidebarMenu(menuItem("Welcome", tabName = "welcome", icon = icon("house")),
                menuItem("Generate database", tabName = "db_generation", icon = icon("database")),
                menuItem("Run proteomics analysis", tabName = "analyse_proteomics", icon = icon("search")),
                menuItem("Integrate data", tabName = "integration", icon = icon("code-merge")),
                menuItem("Visualise results", tabName = "visualisation", icon = icon("eye")),
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

        Shiny.addCustomMessageHandler('switchTab', function(params) {
          var tab = params.tab;
          var target_hash = '#shiny-tab-' + tab;
          var a_tags = document.getElementsByTagName('a');
          for (let i = 0; i < a_tags.length; ++i) {
            let a_tag = a_tags[i];
            let hash = a_tag.hash;
            if (hash === target_hash) {
              a_tag.click();
              break;
            }
          }
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
                           img(src = "images/workflow.png", width = "100%")
                       )
                )
              ),
              fluidRow(
                column(6,
                       h3("Validate installation", style = "margin: 0px"),
                       actionButton("test_data_submit_button", "Run test data", class = "btn btn-info"),
                       downloadButton("test_data_download_button", "Download results (zip)", style = "width:30%;", enabled = FALSE), # initially disabled
                       div(id = "test-data-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "test-data-status-msg-container")
                ),
                column(6,
                       h4("For more information about the test data files used to validate the installation, click the button below to head to the GenomeProt help guide."),
                       actionButton("quick_start_button", "Quick start", class = "btn"),
                       h4("To visualise results from running the test data through GenomeProt, click the button below to head to the visualisation module."),
                       actionButton("visualise_data_button", "Visualise data", class = "btn")
                )
              )
      ),
      tabItem(tabName = "db_generation",
              h2("Generate a custom proteogenomics database"),
              h5("Creates an amino acid FASTA of all ORFs in your data to use as input for FragPipe/MaxQuant etc."),
              fluidRow(
                column(6,
                       # radio buttons for the sequencing type (long-read or short-read) and input type (BAM or GTF)
                       radioButtons("sequencing_type", "Select sequencing type:",
                                    choices = c("Long-read (ONT, PacBio)" = "long-read",
                                                "Short-read" = "short-read")),
                       radioButtons("input_type", "Select input type:",
                                    choices = c("FASTQs" = "fastq_input",
                                                "BAMs" = "bam_input",
                                                "GTF (and/or transcript counts)" = "gtf_input")),

                       # checkbox for whether to generate a variant-aware proteome database
                       checkboxInput("vcf_option", "Incorporate SNVs into protein sequences", value = FALSE),

                       # organism
                       selectInput("organism", "Organism:",
                                   choices = c("Human (H. sapiens)" = "HUMAN",
                                               "Roundworm (C. elegans)" = "CAEEL",
                                               "Fruit fly (D. melanogaster)" = "DROME",
                                               "Mouse (M. musculus)" = "MOUSE",
                                               "Rat (R. norvegicus)" = "RAT",
                                               "Zebrafish (D. rerio)" = "DANRE")),
                                               #"Chimpanzee (P. troglodytes)" = "PANTR",
                                               #"Cow (B. taurus)" = "BOVIN",
                                               #"Clawed frog (X. tropicalis)" = "XENTR",
                                               #"Baker's yeast (S. cerevisiae)" = "YEAST")),

                       # Type of ORFs to include in the proteome database
                       radioButtons("database_type", "Type of ORFs to include in the proteome database:",
                                    choices = c("Canonical" = "canonical",
                                                "All" = "all"),
                                    selected = "all"),

                       # ORF length cutoff
                       numericInput("min_orf_length", "Minimum ORF length (in amino acids):", value = 30, min = 0, step = 1),

                       # options for finding short ORFs
                       h5(tags$b("Find short (10 to 'Minimum ORF length' amino acids) ORFs in the UTRs of reference transcripts:")),
                       checkboxInput("user_find_utr_5_orfs",   "Upstream 5' ORFs"),
                       checkboxInput("user_find_utr_3_orfs", "Downstream 3' ORFs"),
                       numericInput("minimum_tx_count", "Minimum expression threshold (sum per transcript):", value = 5, min = 0, step = 1),
                       numericInput("user_threads",       "CPUs:", value = 1, min = 1, step = 1),

                       # FASTQ-specific input options
                       conditionalPanel(condition = "input.input_type == 'fastq_input'",
                         fileInput("user_reference_genome", "Upload reference genome FASTA file (can be gzipped):",        accept = c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn",
                                                                                                                                      ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")),

                         conditionalPanel(condition = "input.sequencing_type == 'short-read'",
                           fileInput("transcriptome_file",  "Upload reference transcriptome FASTA file (can be gzipped):", accept = c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn",
                                                                                                                                      ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz"))
                         ),

                         fileInput("user_fastq_files",      "Upload FASTQ files (can be gzipped):",                        accept = c(".fastq", ".fq", ".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn",
                                                                                                                                      ".fastq.gz", ".fq.gz", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz"), multiple = TRUE),
                       ),

                       # BAM-specific input options
                       conditionalPanel(condition = "input.input_type == 'bam_input'",
                         fileInput("user_reference_genome", "Upload reference genome FASTA file (can be gzipped):", accept = c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn",
                                                                                                                               ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")),
                         fileInput("user_bam_files",        "Upload BAM files:",                                    accept = c(".bam"), multiple = TRUE),
                       ),

                       # GTF-specific input options
                       conditionalPanel(condition = "input.input_type == 'gtf_input'",
                         conditionalPanel(condition = "input.sequencing_type == 'long-read'",
                           fileInput("user_gtf_file",         "Upload user-generated transcript annotation GTF file (e.g. 'bambu_transcript_annotations.gtf'):", accept = c(".gtf", ".gff", ".gff2", ".gff3")),
                           fileInput("user_tx_count_file",    "Upload user-generated transcript counts file (optional; e.g. 'bambu_transcript_counts.txt'):",    accept = c(".txt", ".csv", ".tsv"))
                         ),

                         conditionalPanel(condition = "input.sequencing_type == 'short-read'",
                           fileInput("user_tx_count_file",    "Upload user-generated transcript counts file (e.g. 'bambu_transcript_counts.txt'):", accept = c(".txt", ".csv", ".tsv"))
                         ),

                         conditionalPanel(condition = "input.vcf_option == true",
                           fileInput("user_reference_genome", "Upload reference genome FASTA file (can be gzipped):", accept = c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn",
                                                                                                                                 ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz"))
                         )
                       ),

                       # reference transcriptome GTF file
                       fileInput("reference_gtf_file", "Upload reference transcriptome GTF file:", accept = c(".gtf", ".gff", ".gff2", ".gff3")),

                       # VCF file
                       conditionalPanel(condition = "input.vcf_option == true",
                         fileInput("user_vcf_file", "Upload VCF file:", accept = c(".vcf"))
                       ),

                       actionButton("db_submit_button", "Submit", class = "btn btn-info")
                ),
                column(6,
                       h3("Download your results:"),
                       downloadButton("db_download_button", "Download results (zip)", style = "width:70%;", enabled = FALSE), # initially disabled
                       div(id = "db-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "db-status-msg-container")
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
                       fileInput("fragpipe_prot_db_fasta_file", "Upload a proteome database FASTA file (e.g. 'proteome_database.fasta'):", accept = c(".fasta")),
                       checkboxInput("user_add_contaminants", "Add contaminants into the proteome database?", value = TRUE),
                       checkboxInput("user_perform_quantification", "Perform peptide quantification after the peptide search?"),
                       h5(tags$b("Upload mass spectrometry data files and select their data types:")),
                       actionButton("add_mass_spec_file_button", "+ Add mass spectrometry data file", class = "btn btn-warning"),
                       actionButton("remove_mass_spec_file_button", "- Remove mass spectrometry data file", class = "btn"),
                       div(id = "mass_spec_file_list"),
                       selectInput("protease1", "Protease 1:",
                                   choices = c("stricttrypsin (cuts KR, sense C)" = "stricttrypsin",
                                               "trypsin (cuts KR, no cuts P, sense C)" = "trypsin",
                                               "trypsin_gluc (cuts DEKR, no cuts P, sense C)" = "trypsin_gluc",
                                               "gluc (cuts DE, no cuts P, sense C)" = "gluc",
                                               "lysc (cuts K, no cuts P, sense C)" = "lysc",
                                               "lysn (cuts K, sense N)" = "lysn",
                                               "argc (cuts R, no cuts P, sense C)" = "argc",
                                               "aspn (cuts D, sense N)" = "aspn")),
                       selectInput("protease2", "Protease 2 (optional; must be different from protease 1):",
                                   choices = c("None" = "none",
                                               "stricttrypsin (cuts KR, sense C)" = "stricttrypsin",
                                               "trypsin (cuts KR, no cuts P, sense C)" = "trypsin",
                                               "trypsin_gluc (cuts DEKR, no cuts P, sense C)" = "trypsin_gluc",
                                               "gluc (cuts DE, no cuts P, sense C)" = "gluc",
                                               "lysc (cuts K, no cuts P, sense C)" = "lysc",
                                               "lysn (cuts K, sense N)" = "lysn",
                                               "argc (cuts R, no cuts P, sense C)" = "argc",
                                               "aspn (cuts D, sense N)" = "aspn")),
                       numericInput("fragpipe_cpu_threads", "Number of CPU threads to use (specify '0' for FragPipe to use <number of cores in system - 1> threads)", value = 1, min = 0, step = 1),
                       numericInput("fragpipe_memory_limit", "Memory limit in GB to use (specify '0' to let FragPipe decide)", value = 15, min = 0, step = 1),
                       actionButton("fragpipe_submit_button", "Run FragPipe", class = "btn btn-info")
                ),
                column(6,
                       h3("Download FragPipe results:"),
                       downloadButton("fragpipe_download_button", "Download FragPipe results (zip)", style = "width:70%;", enabled = FALSE), # initially disabled
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
                       h3("Part 1: Reformat proteomics results files"),
                       h5("Note 1: All proteomics results files with a file extension of '.txt' or '.csv' will be renamed to have a flie extension of '.tsv'."),
                       h5("Note 2: Ignoring the file extension, if a proteomics results file with the name 'peptide_data' was uploaded, it will be renamed to 'peptide_data_renamed.tsv'."),
                       fileInput("user_orig_proteomics_files", "Upload proteomics results files:", accept = c(".txt", ".csv", ".tsv"), multiple = TRUE),
                       radioButtons("proteomics_search_tool", "Select proteomics search tool:",
                                    choices = c("Spectronaut" = "Spectronaut",
                                                "FragPipe (identified peptides, i.e. peptide.tsv)" = "FragPipe",
                                                "FragPipe (quantified peptides, i.e. report.pr_matrix.tsv)" = "FragPipe_quant")),
                       actionButton("integ_reformat_submit_button", "Submit", class = "btn btn-info")
                ),
                column(6,
                       h3("Download the reformatted proteomics results file:"),
                       downloadButton("integ_reformat_download_button", "Download reformatted results file (peptide_data.tsv)", style = "width:70%;", enabled = FALSE), # initially disabled
                       div(id = "integ-reformat-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "integ-reformat-status-msg-container")
                )
              ),
              fluidRow(
                column(6,
                       h3("Part 2: Upload files to integrate"),
                       fileInput("user_proteomics_file", "Upload reformatted proteomics results file:", accept = c(".txt", ".csv", ".tsv")),
                       fileInput("user_metadata_file",   "Upload 'proteome_database_metadata.txt':",    accept = c(".txt")),
                       fileInput("user_post_gtf_file",   "Upload 'proteome_database_transcripts.gtf':", accept = c(".gtf")),
                       actionButton("integ_submit_button", "Submit", class = "btn btn-info")
                ),
                column(6,
                       h3("Download integration results:"),
                       downloadButton("integ_download_button", "Download results (zip)", style = "width:70%;", enabled = FALSE), # initially disabled
                       div(id = "integ-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "integ-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "visualisation",
              fluidRow(
                column(12,
                    tags$iframe(id = "isovis_window",
                                src = "https://isomix.org/isovis/",
                                width = "100%",
                                height = "950px",
                                style = "border:none;"))
              ),
              h5("The IsoVis website is displayed above for convenience. It is also accessible directly at: https://isomix.org/isovis/"),
              h5(actionLink("show_isovis_steps", "Instructions for using IsoVis")),
              conditionalPanel(
                condition = "input.show_isovis_steps % 2 == 1",
                p("Step 1: Click 'Upload data'."),
                p("Step 2: For the 'Transcript data' file, upload 'combined_annotations.gtf'. For the 'Transcript counts' file, upload 'bambu_transcript_counts.txt'."),
                p("Step 3: For the 'Peptide intensities' file, upload the peptide intensities file from the proteomics pipeline you used (e.g. 'report.pr_matrix.tsv'), then click 'Apply'."),
                p("Step 4: Type the symbol or ID of a gene to view, select it from the list of results displayed, then either press enter or click '>'."),
                p("Step 5: To see the mappings of peptides to open reading frames, click on the 'Stack options' dropdown menu and select 'Peptide mapping'.")
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