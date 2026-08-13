library(shiny)
library(shinyjs)
library(reticulate)

conda_path <-conda_binary()
conda_command <- paste0(conda_path, " run -n GenomeProt_env ")

# TODO: R_ZIPCMD might be defined but empty sometimes, which causes zip() to fail
if (nchar(Sys.getenv("R_ZIPCMD", "zip")) == 0) {
  Sys.setenv("R_ZIPCMD" = "zip")
}

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
    command_gffread <- paste0(conda_command, "gffread -w ", transcript_fasta_file, " -g ", genome_fasta_file, " ", input$reference_gtf_file$datapath)
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
    logfile_path <- file.path(outdir_bambu, "logfile.txt")

    if ( !is.null(input$bam_data) && input$bam_data == "user" && !is.null(input$user_bam_files)) {
      bam_df <- as.data.frame(input$user_bam_files)
      bamdir <- dirname(input$user_bam_files$datapath)

      # Specify new filenames
      new_names <- file.path(bamdir[1], input$user_bam_files$name)

      # Rename the files
      file.rename(input$user_bam_files$datapath, new_names)

      command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R -b ", bamdir[1],
                              " -g ", input$reference_gtf_file$datapath, " -o ", outdir_bambu, " -t ", input$user_threads, " -s ", input$organism, " > ", logfile_path, " 2>&1")
    } else {
      command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R -b ", default_bam_path,
                              " -g ",                  default_gtf_path, " -o ", outdir_bambu, " -t ", input$user_threads, " -s ", input$organism, " > ", logfile_path, " 2>&1")
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
    req(input$user_gtf_file) # GTFs required
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
    ref_genome <- "./testdata/GRCh38_chr1_6_7_masked.fa.gz"
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
    ref_proteome <- "data/openprot_uniprotDb_human.txt"
  } else if (input$organism == "MOUSE") {
    ref_proteome <- "data/openprot_uniprotDb_mouse.txt"
  } else if (input$organism == "CAEEL") {
    ref_proteome <- "data/openprot_uniprotDb_c_elegans.txt"
  } else if (input$organism == "DROME") {
    ref_proteome <- "data/openprot_uniprotDb_drosophila.txt"
  } else if (input$organism == "RAT") {
    ref_proteome <- "data/openprot_uniprotDb_rat.txt"
  } else if (input$organism == "DANRE") {
    ref_proteome <- "data/openprot_uniprotDb_zebrafish.txt"
  }else if (input$organism == "PANTR") {
    ref_proteome <- "data/openprot_uniprotDb_chimp.txt"
  }else if (input$organism == "BOVIN") {
    ref_proteome <- "data/openprot_uniprotDb_cow.txt"
  }else if (input$organism == "XENTR") {
    ref_proteome <- "data/openprot_uniprotDb_clawed_frog.txt"
  }else if (input$organism == "YEAST") {
    ref_proteome <- "data/openprot_uniprotDb_yeast.txt"
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
      files_to_zip_db <- c("../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/bambu_transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "../bambu_output/logfile.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
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

# Proteomics module
is_fragpipe_input_valid <- function(input, mass_spec_file_num) {
  error_msg <- ""

  # Check if the user has uploaded a proteome database
  fragpipe_prot_db_file <- input$fragpipe_prot_db_fasta_file
  if (is.null(fragpipe_prot_db_file)) {
    error_msg <- "Error: Please upload a proteome database FASTA file to perform a proteomics search."
    return(error_msg)
  }

  # The proteome database must have a file extension of '.fasta'
  filename <- fragpipe_prot_db_file$name
  if (!endsWith(filename, ".fasta")) {
    error_msg <- "Error: The proteome database file must have a file extension of '.fasta'."
    return(error_msg)
  }

  # Check if the user has uploaded at least one mass spectrometry data file
  if (mass_spec_file_num == 1) {
    error_msg <- "Error: Please upload at least one mass spectrometry data file to perform a proteomics search."
    return(error_msg)
  }

  # Check if the user uploaded exactly (mass_spec_file_num - 1) mass spectrometry data files and selected that number of data types
  for (i in 1:(mass_spec_file_num - 1)) {
    mass_spec_file <- input[[paste0("mass_spec_file_input_", i)]]
    if (is.null(mass_spec_file)) {
      error_msg <- "Error: Please ensure all mass spectrometry data file inputs have an uploaded file."
      return(error_msg)
    }

    # Mass spectrometry data files must have a file extension of '.mzML', '.mzXML', '.mgf', '.mzBIN', '.raw', or '.d'
    filename <- mass_spec_file$name
    if (!base::any(endsWith(filename, c(".mzML", ".mzXML", ".mgf", ".mzBIN", ".raw", ".d")))) {
      error_msg <- "Error: All mass spectrometry data files must have a file extension of '.mzML', '.mzXML', '.mgf', '.mzBIN', '.raw' or '.d'."
      return(error_msg)
    }

    # Mass spectrometry data files must have a data type of 'DDA', 'DDA+', 'DIA', 'DIA-Quant', 'DIA-Lib', or 'GPF-DIA'
    mass_spec_data_type <- input[[paste0("mass_spec_file_data_type_input_", i)]]
    if (is.null(mass_spec_data_type) || !(mass_spec_data_type %in% c("DDA", "DDA+", "DIA", "DIA-Quant", "DIA-Lib", "GPF-DIA"))) {
      error_msg <- "Error: All mass spectrometry data files must have a data type of 'DDA', 'DDA+', 'DIA', 'DIA-Quant', 'DIA-Lib' or 'GPF-DIA'."
      return(error_msg)
    }
  }

  # The first protease used must be 'stricttrypsin', 'trypsin', 'trypsin_gluc', 'gluc', 'lysc', 'lysn', 'argc', or 'aspn'
  protease1 <- input$protease1
  if (is.null(protease1) || !(protease1 %in% c("stricttrypsin", "trypsin", "trypsin_gluc", "gluc", "lysc", "lysn", "argc", "aspn"))) {
    error_msg <- "Error: The first protease used must be 'stricttrypsin', 'trypsin', 'trypsin_gluc', 'gluc', 'lysc', 'lysn', 'argc', or 'aspn'."
    return(error_msg)
  }

  # If the second protease were specified, it must be 'stricttrypsin', 'trypsin', 'trypsin_gluc', 'gluc', 'lysc', 'lysn', 'argc', or 'aspn'
  protease2 <- input$protease2
  if (!(is.null(protease2) || (protease2 == "none")) && !(protease2 %in% c("stricttrypsin", "trypsin", "trypsin_gluc", "gluc", "lysc", "lysn", "argc", "aspn"))) {
    error_msg <- "Error: The second protease used must be 'stricttrypsin', 'trypsin', 'trypsin_gluc', 'gluc', 'lysc', 'lysn', 'argc', or 'aspn'."
    return(error_msg)
  }

  # The first protease must be different from the second protease
  if (protease1 == protease2) {
    error_msg <- "Error: The first protease must be different from the second protease."
    return(error_msg)
  }

  # Number of CPU threads: Integer >= 0
  fragpipe_cpu_threads <- floor(input$fragpipe_cpu_threads)
  if (!(is.finite(fragpipe_cpu_threads) && (fragpipe_cpu_threads >= 0))) {
    error_msg <- paste0("Error: The number of CPU threads must be a non-negative integer. Entered value: ", fragpipe_cpu_threads)
    return(error_msg)
  }

  # Memory limit: Integer >= 0
  fragpipe_memory_limit <- floor(input$fragpipe_memory_limit)
  if (!(is.finite(fragpipe_memory_limit) && (fragpipe_memory_limit >= 0))) {
    error_msg <- paste0("Error: The memory limit in GB must be a non-negative integer. Entered value: ", fragpipe_memory_limit)
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

fragpipe_server <- function(input, output, session_id, mass_spec_file_num) {
  fragpipe_prot_db_file <- input$fragpipe_prot_db_fasta_file
  renamed_fragpipe_prot_db_file <- file.path(dirname(fragpipe_prot_db_file$datapath), fragpipe_prot_db_file$name)
  if (file.exists(fragpipe_prot_db_file$datapath)) {
    file.rename(from = fragpipe_prot_db_file$datapath, to = renamed_fragpipe_prot_db_file)
  }

  mass_spec_info_file_contents <- ''
  for (i in 1:(mass_spec_file_num - 1)) {
    mass_spec_file <- input[[paste0("mass_spec_file_input_", i)]]
    mass_spec_data_type <- input[[paste0("mass_spec_file_data_type_input_", i)]]

    renamed_file <- file.path(dirname(mass_spec_file$datapath), mass_spec_file$name)
    if (file.exists(mass_spec_file$datapath)) {
      file.rename(from = mass_spec_file$datapath, to = renamed_file)
    }

    mass_spec_info_file_contents <- paste0(mass_spec_info_file_contents, renamed_file, "\t", mass_spec_data_type, "\n")
  }

  # Remove the last newline from the mass spectrometry file list contents
  mass_spec_info_file_contents <- substr(mass_spec_info_file_contents, 1, nchar(mass_spec_info_file_contents) - 1)

  # Go into the session directory
  old_wd <- getwd()
  setwd(session_id)

  protease1 <- input$protease1
  protease2 <- input$protease2
  if (is.null(protease2)) {
    protease2 <- "none"
  }

  fragpipe_cpu_threads <- floor(input$fragpipe_cpu_threads)
  fragpipe_memory_limit <- floor(input$fragpipe_memory_limit)
  output_dir <- "fragpipe_output"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Go into the output directory so that the FragPipe log file will be output there
  output_dir <- normalizePath(output_dir)
  setwd(output_dir)

  # Write the mass spectrometry file list
  mass_spec_info_file <- "mass_spec_info_list.txt"
  writeLines(mass_spec_info_file_contents, mass_spec_info_file)

  command_run_fragpipe <- paste0(
    conda_command, "--no-capture-output python ../../bin/proteomics_module/fragpipe-run.py",
    " --db_path ", shQuote(renamed_fragpipe_prot_db_file),
    " --mass_spec_info_path ", shQuote(mass_spec_info_file),
    " --output_dir ", shQuote(output_dir),
    " --protease1 ", shQuote(protease1)
  )

  if (protease2 != "none") {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --protease2 ", shQuote(protease2))
  }

  if (input$user_add_contaminants) {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --add_contaminants")
  }

  if (input$user_perform_quantification) {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --perform_quantification")
  }

  command_run_fragpipe <- paste0(
    command_run_fragpipe,
    " --fragpipe_path ", shQuote("/home/user/Desktop/GenomeProt/fragpipe-23.1/"),
    " --num_threads ", shQuote(fragpipe_cpu_threads),
    " --memory_limit ", shQuote(fragpipe_memory_limit)
  )
  message(command_run_fragpipe)
  system(command_run_fragpipe)

  # Find files to zip
  error_msg <- ''

  files_to_zip_fragpipe <- c("peptide.tsv")
  if (input$user_perform_quantification) {
    files_to_zip_fragpipe <- c(files_to_zip_fragpipe, file.path("dia-quant-output", "report.pr_matrix.tsv"))
  }

  if (base::all(file.exists(files_to_zip_fragpipe))) {
    zipfile_path_fragpipe <- file.path("..", "fragpipe_results.zip")
    zip(zipfile = zipfile_path_fragpipe, files = files_to_zip_fragpipe)
  } else {
    error_msg <- "MISSING_RESULTS"
  }

  # Go back to the old directory
  setwd(old_wd)

  return(error_msg)
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

  # PROTEOMICS MODULE

  # create reactive value for the FragPipe results zip file
  file_available_fragpipe <- reactiveVal(FALSE)

  # add a counter for the mass spectrometry files
  mass_spec_file_num <- 1

  hide("remove_mass_spec_file_button")

  # add a mass spectrometry file upload button and its corresponding data type dropdown menu when
  # the '+ Add mass spectrometry data file' button is pressed
  observeEvent(input$add_mass_spec_file_button, {
    insertUI(
      selector = "#mass_spec_file_list",
      ui = fluidRow(id = paste0("mass_spec_row_", mass_spec_file_num),
              column(8,
                     fileInput(paste0("mass_spec_file_input_", mass_spec_file_num),
                               label = paste0("Upload mass spectrometry file #", mass_spec_file_num, ":"),
                               buttonLabel = "Browse...", multiple = FALSE, accept = c(".mzML", ".mzXML", ".mgf", ".mzBIN", ".raw", ".d"))
              ),
              column(4,
                     selectInput(paste0("mass_spec_file_data_type_input_", mass_spec_file_num),
                                 label = "Data type:",
                                 choices = list("DDA" = "DDA",
                                                "DDA+" = "DDA+",
                                                "DIA" = "DIA",
                                                "DIA-Quant" = "DIA-Quant",
                                                "DIA-Lib" = "DIA-Lib",
                                                "GPF-DIA" = "GPF-DIA"),
                                 selected = "DDA")
              )
           )
    )
    shinyjs::runjs("document.getElementById('mass_spec_file_list').style = 'border: 2px solid; padding: 5px;';")
    shinyjs::runjs("document.getElementById('add_mass_spec_file_button').children[0].innerText = '+ Add another mass spectrometry data file'")
    mass_spec_file_num <<- mass_spec_file_num + 1
    showElement("remove_mass_spec_file_button")
  })

  # add a mass spectrometry file upload button and its corresponding data type dropdown menu when
  # the '+ Add mass spectrometry data file' button is pressed
  observeEvent(input$remove_mass_spec_file_button, {
    mass_spec_file_num <<- mass_spec_file_num - 1
    removeUI(selector = paste0("#mass_spec_row_", mass_spec_file_num))
    if (mass_spec_file_num == 1) {
      shinyjs::runjs("document.getElementById('mass_spec_file_list').style = '';")
      shinyjs::runjs("document.getElementById('add_mass_spec_file_button').children[0].innerText = '+ Add mass spectrometry data file'")
      hide("remove_mass_spec_file_button")
    }
  })

  # run FragPipe when the 'Run FragPipe' button is pressed
  observeEvent(input$fragpipe_submit_button, {
    file_available_fragpipe(FALSE)

    message("Validating FragPipe input...")
    error_msg <- is_fragpipe_input_valid(input, mass_spec_file_num)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "fragpipe-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "fragpipe-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "fragpipe_submit_button", spinnerId = "fragpipe-loading-container")) # disable submit button

    # run the FragPipe server
    if (file.exists(file.path(session_id, "fragpipe_results.zip"))) {
      file.remove(file.path(session_id, "fragpipe_results.zip"))
    }
    error_msg <- fragpipe_server(input, output, session_id, mass_spec_file_num)

    session$sendCustomMessage("enableButton", list(id = "fragpipe_submit_button", spinnerId = "fragpipe-loading-container"))

    # check if the zip file is created
    if ((error_msg != '') || (!file.exists(file.path(session_id, "fragpipe_results.zip")))) {
      error_msg <- "Error: Proteomics search failed. Please check the proteomics module log file generated in the output directory for details."
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "fragpipe-status-msg-container", color = "red"))
    } else {
      file_available_fragpipe(TRUE)
    }
  })

  # enable download once files are available
  observe({
    if (file_available_fragpipe()) {
      shinyjs::enable("fragpipe_download_button")
      shinyjs::runjs("document.getElementById('fragpipe_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "fragpipe_submit_button", spinnerId = "fragpipe-loading-container")) # re-enable submit button
    }
  })

  # download handler for the FragPipe results zip
  output$fragpipe_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_fragpipe_results.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "fragpipe_results.zip"), file)
    }
  )

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