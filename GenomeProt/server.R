library(shiny)
library(shinyjs)
library(reticulate)

conda_path <- conda_binary()
conda_command <- paste0(conda_path, " run -n GenomeProt_env ")

# internal server functions
bam_server <- function(input, output, session) {
  req(input$user_reference_genome_bam$datapath,input$reference_gtf_file$datapath)   # BAMs and ref GTF required

  session_id <- session$token   # session ID
  outdir_bam <- file.path(session_id, "mapping_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_bam)) {
    dir.create(outdir_bam)
  }

  # generate reference.fa
  genome_fasta_file     <- file.path(outdir_bam,     "genome.fa")
  transcript_fasta_file <- file.path(outdir_bam, "transcript.fa")
  if (grepl("\\.gz$", input$user_reference_genome_bam$datapath)) {  # check if file is compressed
    command_decompress <- paste0(conda_command, "gzip -c -d ", input$user_reference_genome_bam$datapath, " >", genome_fasta_file)
    print(command_decompress)
    system(command_decompress)
    command_gffread <- paste0(conda_command, "gffread -w ", transcript_fasta_file, " -g ",                        genome_fasta_file, " ", input$reference_gtf_file$datapath)
  } else {
    command_gffread <- paste0(conda_command, "gffread -w ", transcript_fasta_file, " -g ", input$user_reference_genome_bam$datapath, " ", input$reference_gtf_file$datapath)
  }

  print(command_gffread)
  system(command_gffread)

  # create df of bam file names
  user_bam_files_df <- input$user_bam_files %>%
    mutate(file_prefix = sub("\\.bam$", "", name))

  # for each bam file
  for (i in 1:nrow(user_bam_files_df)) {
    bam_file <- user_bam_files_df$datapath[i]
    file_prefix <- user_bam_files_df$file_prefix[i]

    command_salmon <- paste0(conda_command, "salmon quant -t ", transcript_fasta_file, " -p ", input$user_threads, " -l A -a ", bam_file, " -o ", file.path(outdir_bam, file_prefix))

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
  gtf_data <- rtracklayer::import(input$reference_gtf_file$datapath, format = "gtf")

  # convert GTF data to a data frame
  gtf_df <- as.data.frame(gtf_data)

  # filter for relevant columns
  transcript_gene_info <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
  colnames(transcript_gene_info) <- c("TXNAME", "GENEID")

  # convert counts object to df
  count_df <- as.data.frame(txi$counts)

  # set col names
  count_df <- count_df %>% mutate(TXNAME = rownames(count_df)) %>% dplyr::select(TXNAME, everything())
  rownames(count_df) <- NULL

  # merge tx counts and gene info
  count_df_merged <- left_join(count_df, transcript_gene_info, by = "TXNAME")
  count_df_merged <- count_df_merged %>% dplyr::select(TXNAME, GENEID, everything())

  # export short-read count data
  write_tsv(count_df_merged, file = file.path(outdir_bam, "bambu_transcript_counts.txt"), escape = "none", col_names = TRUE)
}

bambu_server <- function(input, output, session) {
  session_id <- session$token   # session ID
  outdir_bambu <- file.path(session_id, "bambu_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_bambu)) {
    dir.create(outdir_bambu)
  }

  if (input$input_type == "bam_input") { # if user input bam files
    # if provided, use user-uploaded files instead of the default files
    default_gtf_path <- "./testdata/gencode_v47_sorted.gtf"
    default_bam_path <- "./testdata/long_read_bam/"

    if ( !is.null(input$bam_data) && input$bam_data == "user" && !is.null(input$user_bam_files)) {
      bam_df <- as.data.frame(input$user_bam_files)
      bamdir <- dirname(input$user_bam_files$datapath)

      # Specify new filenames
      new_names <- file.path(bamdir[1], input$user_bam_files$name)

      # Rename the files
      file.rename(input$user_bam_files$datapath, new_names)

      command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R -b ", bamdir[1],
                             "/ -g ", input$reference_gtf_file$datapath, " -o ", outdir_bambu, " -t ", input$user_threads, " -s ", input$organism)
    } else {
      command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R -b ", default_bam_path,
                             "/ -g ",                  default_gtf_path, " -o ", outdir_bambu, " -t ", input$user_threads, " -s ", input$organism)
    }
  }

  # run bambu
  print(command_bambu)
  system(command_bambu)

  # rename bambu output files
  renamed_gtf <- file.path(outdir_bambu, "bambu_transcript_annotations.gtf")
  file.rename(file.path(outdir_bambu,    "counts_transcript.txt"), file.path(outdir_bambu, "bambu_transcript_counts.txt"))
  file.rename(file.path(outdir_bambu, "extended_annotations.gtf"), renamed_gtf)

  if (!is.null(input$bam_data) && input$bam_data == "user" && !is.null(input$user_bam_files)) {
    command_gff_compare <- paste0(conda_command, "gffcompare -r ", input$reference_gtf_file$datapath, " ", renamed_gtf)
  } else {
    command_gff_compare <- paste0(conda_command, "gffcompare -r ",                  default_gtf_path, " ", renamed_gtf)
  }

  # run gffcompare
  print(command_gff_compare)
  system(command_gff_compare)

  # rename gffcompare output files
  file.rename(file.path(outdir_bambu, "gffcmp.bambu_transcript_annotations.gtf.tmap"), file.path(outdir_bambu, "gffcompare.tmap.txt"))
  file.remove(Sys.glob("gffcmp*"))
}

database_server <- function(input, output, session) {
  session_id <- session$token   # session ID
  outdir_db <- file.path(session_id, "database_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_db)) {
    dir.create(outdir_db)
  }

  # set input file type
  if (input$input_type == "gtf_input") {
    req(input$user_gtf_file, input$reference_gtf_file) # GTFs required
    db_gtf_file <- input$user_gtf_file$datapath
    db_counts_file <- input$user_tx_count_file$datapath
  } else if (input$input_type == "bam_input") {
    if (input$sequencing_type == "long-read") {
      db_gtf_file <- file.path(session_id, "bambu_output", "bambu_transcript_annotations.gtf")
      db_counts_file <- file.path(session_id, "bambu_output", "bambu_transcript_counts.txt")
    } else if (input$sequencing_type == "short-read") {
      db_gtf_file <- input$reference_gtf_file$datapath
      db_counts_file <- file.path(session_id, "mapping_output", "bambu_transcript_counts.txt")
    }
  }

  if (!is.null(input$reference_gtf_file$datapath) && nzchar(input$reference_gtf_file$datapath)) {
    ref_gtf <- input$reference_gtf_file$datapath
  } else {
    ref_gtf <- "./testdata/gencode_v47_sorted.gtf"
  }

  # conditionally add counts file argument
  counts_arg <- if (!is.null(db_counts_file) && nzchar(db_counts_file)) {
    paste0(" -c ", db_counts_file)
  } else {
    ""
  }

  # conditionally add genome input argument
  genome_arg <- if (!is.null(input$user_reference_genome_bam$datapath)) {
    paste0(" -G ", input$user_reference_genome_bam$datapath)
  } else {
    ref_genome <- "./testdata/GRCh38_chr1_6_7.fa"
    paste0(" -G ", ref_genome)
  }
  
  vcf_arg <- if (input$vcf_option == TRUE && !is.null(input$user_vcf_file) && !is.null(genome_arg)) {
    vcf_file <- input$user_vcf_file
    paste0(" -v ", vcf_file, genome_arg)
  } else if (input$vcf_option == TRUE && is.null(input$user_vcf_file)) {
    vcf_file <- "./testdata/BRAF_mutation.vcf"
    paste0(" -v ", vcf_file, genome_arg)
  } else {
    vcf_file <- NULL
    ""
  }

  # construct the command
  command_generate_proteome <- paste0(
    conda_command, "Rscript bin/database_module/generate_proteome.R",
    " -g ", db_gtf_file,
    " -r ", ref_gtf,
    counts_arg, # include counts file only if provided
    " -m ", input$minimum_tx_count,
    " -o ", input$organism,
    " -l ", input$min_orf_length,
    " -u ", input$user_find_utr_5_orfs,
    " -d ", input$user_find_utr_3_orfs,
    vcf_arg,    # include VCF file only if provided
    " -s ", outdir_db
  )

  print(genome_arg)
  print(input$vcf_option)
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
  orfome_file <- file.path(outdir_db, "ORFome_aa.txt")
  proteome_database_transcripts_gtf_file <- file.path(outdir_db, "proteome_database_transcripts.gtf")
  mutant_orfome_file <- file.path(outdir_db, "Mutant_ORFome_aa.txt")
  if (!is.null(vcf_file)) { # if there is a VCF file uploaded
    command_annotate_proteome <- paste0(conda_command, "--no-capture-output python bin/database_module/annotate_proteome.py ", ref_gtf, " ", ref_proteome, " ", orfome_file, " ", proteome_database_transcripts_gtf_file, " ",
                                        outdir_db, " ", input$database_type, " ", input$min_orf_length, " ", mutant_orfome_file, " ", input$organism, " ", input$user_threads, " 2000")
  } else { # if no VCF file uploaded
    command_annotate_proteome <- paste0(conda_command, "--no-capture-output python bin/database_module/annotate_proteome.py ", ref_gtf, " ", ref_proteome, " ", orfome_file, " ", proteome_database_transcripts_gtf_file, " ",
                                        outdir_db, " ", input$database_type, " ", input$min_orf_length, " ",             "None", " ", input$organism, " ", input$user_threads, " 2000")
  }

  # run python script to create proteome fasta
  print(command_annotate_proteome)
  system(command_annotate_proteome)
  print("Annotated proteome")

  # get top level directory
  top_level_dir <- getwd()

  # zip all results files depending on input types
  proteome_database_fasta_file <- file.path(outdir_db, "proteome_database.fasta")
  orf_temp_file <- file.path(outdir_db, "orf_temp.txt")
  if (file.exists(proteome_database_fasta_file) && file.exists(proteome_database_transcripts_gtf_file) && !file.exists(orf_temp_file)) {
    if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      files_to_zip_db <- c("../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/bambu_transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    } else if (input$sequencing_type == "short-read") {
      files_to_zip_db <- c("../mapping_output/bambu_transcript_counts.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    } else if (input$input_type == "gtf_input") {
      files_to_zip_db <- c("proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    }

    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_db <- file.path("..", "database_results.zip")

    # temp change the working dir to outdir_db
    tmp_wd <- setwd(outdir_db)

    if (file.exists(zipfile_path_db)) {
      file.remove(zipfile_path_db)
    }

    # zip files
    zip(zipfile = zipfile_path_db, files = files_to_zip_db)

    # change back to starting wd
    setwd(top_level_dir)
  }
}

integration_server <- function(input, output, session) {
  req(input$user_proteomics_file, input$user_post_gtf_file, input$user_metadata_file)  # GTF is required

  session_id <- session$token   # session ID
  outdir_integ <- file.path(session_id, "integ_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_integ)) {
    dir.create(outdir_integ)
  }

  # run Rscript
  system(paste0(conda_command, "Rscript bin/integration_module/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -m ", input$user_metadata_file$datapath, " -g ", input$user_post_gtf_file$datapath, " -s ", outdir_integ))

  # get the top level dir
  top_level_dir <- getwd()

  # create report
  report_rmd_file <- file.path(top_level_dir, "bin", "integration_module", "integration_summary_report.Rmd")
  report_outdir <- file.path(top_level_dir, outdir_integ)
  system(paste0(conda_command, "Rscript bin/integration_module/generate_integration_summary_report.R -i ", report_rmd_file, " -o ", report_outdir))

  # zip all results files
  if (file.exists(file.path(outdir_integ, "peptide_info.tsv")) && file.exists(file.path(outdir_integ, "summary_report.html"))) {
    # create a zip file with results
    files_to_zip_int <- c("summary_report.html", "peptide_info.tsv", "report_images/",
                          "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                          "peptides.bed12", "ORFs.bed12", "transcripts.bed12", "ncORF_stats.xlsx")

    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_int <- file.path("..", "integration_results.zip")

    # temp change the working dir to outdir_integ
    tmp_wd <- setwd(outdir_integ)

    if (file.exists(zipfile_path_int)) {
      file.remove(zipfile_path_int)
    }

    # zip files
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)

    # go back to starting dir
    setwd(top_level_dir)
  }
}

# main shiny app server
server <- function(input, output, session) {
  session_id <- session$token   # session ID
  print(paste0("Session: ", session_id))

  # create session id tmp directory each time app is run
  if (!dir.exists(session_id)) {
    dir.create(session_id)
  }

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
    if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      bambu_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "bam_input" & input$sequencing_type == "short-read") {
      bam_server(input, output, session)
      database_server(input, output, session)
    } else if (input$input_type == "gtf_input") {
      database_server(input, output, session)
    }

    # check if the zip file is created
    if (file.exists(file.path(session_id, "database_results.zip"))) {
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
      file.copy(file.path(session_id, "database_results.zip"), file)
    }
  )

  # END DATABASE MODULE

  # INTEGRATION MODULE

  # create reactive value for the database zip
  file_available_integ <- reactiveVal(FALSE)

  # run integration function when submit is pressed
  observeEvent(input$integ_submit_button, {
    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button

    # run integration server
    integration_server(input, output, session)

    # check if the zip file is created
    if (file.exists(file.path(session_id, "integration_results.zip"))) {
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
      file.copy(file.path(session_id, "integration_results.zip"), file)
    }
  )

  # END INTEGRATION MODULE

  # remove session id tmp directory created each time app is run
  session$onSessionEnded(function() {
    if (dir.exists(session_id)) {
      unlink(session_id, recursive = TRUE)
    }
  })
}