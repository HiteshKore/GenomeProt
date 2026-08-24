conda_path <- reticulate::conda_binary()
conda_command <- paste0(conda_path, " run -n GenomeProt_env ")

# TODO: R_ZIPCMD might be defined but empty sometimes, which causes zip() to fail
if (nchar(Sys.getenv("R_ZIPCMD", "zip")) == 0) {
  Sys.setenv("R_ZIPCMD" = "zip")
}

# internal server functions
is_test_data_input_valid <- function() {
  error_msg <- ""

  # Check if the test data files are unzipped
  test_data_files <- file.path("testdata",
                               c("BRAF_mutation.vcf", "gencode_v47_sorted.gtf", "GRCh38_chr1_6_7_masked.fa", file.path("long_read_bam", "Melanoma_data_subset.bam"), "peptide_data.tsv"))
  if (!base::all(file.exists(test_data_files))) {
    error_msg <- "Error: Not all test data files are unzipped. Please unzip the test data files in the testdata/ directory by running 'bash prepare_ref_and_test_data.sh' in the terminal."
    return(error_msg)
  }

  # Check if the reference data files are unzipped
  reference_data_files <- file.path("data",
                                    c("openprot_uniprotDb_human.txt", "openprot_uniprotDb_c_elegans.txt", "openprot_uniprotDb_drosophila.txt", "openprot_uniprotDb_mouse.txt", "openprot_uniprotDb_rat.txt", "openprot_uniprotDb_zebrafish.txt"))
  if (!base::all(file.exists(reference_data_files))) {
    error_msg <- "Error: Not all reference data files are unzipped. Please unzip the reference data files in the data/ directory by running 'bash prepare_ref_and_test_data.sh' in the terminal."
    return(error_msg)
  }

  # Check if the GenomeProt_env conda environment is set up
  command_conda_env_list <- paste0(conda_path, " env list")
  command_output <- system(command_conda_env_list, intern = TRUE)
  if (!any(grepl("GenomeProt_env", command_output, fixed = TRUE))) {
    error_msg <- "Error: The GenomeProt_env conda environment is not set up. Please build the conda environment by running either 'conda env create -f conda_env.yaml' (if you wish to use the proteomics module) or 'conda env create -f conda_env_no_fragpipe.yaml' (if you do not wish to use the proteomics module) in the terminal."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

test_data_server <- function(session) {
  session_id <- session$token                                   # session ID
  outdir_test_data <- file.path(session_id, "test_data_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_test_data)) {
    dir.create(outdir_test_data)
  }

  error_msg <- ""

  # input files
  bam_file <- file.path("testdata", "long_read_bam", "Melanoma_data_subset.bam")
  reference_genome_file <- file.path("testdata", "GRCh38_chr1_6_7_masked.fa")
  reference_transcriptome_file <- file.path("testdata", "gencode_v47_sorted.gtf")
  vcf_file <- file.path("testdata", "BRAF_mutation.vcf")
  reformatted_peptide_results_file <- file.path("testdata", "peptide_data.tsv")

  ### bambu server ###
  outdir_test_data_bambu <- file.path(outdir_test_data, "bambu_output")
  if (!dir.exists(outdir_test_data_bambu)) {
    dir.create(outdir_test_data_bambu)
  }

  logfile_path <- file.path(outdir_test_data_bambu, "logfile.txt")
  command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R",
                          " -b ", shQuote(dirname(bam_file)),
                          " -g ", shQuote(reference_transcriptome_file),
                          " -o ", shQuote(outdir_test_data_bambu),
                          " -t ", shQuote(1),
                          " -s ", shQuote("HUMAN"),
                          " > ", shQuote(logfile_path),
                          " 2>&1")

  # run bambu
  message(command_bambu)
  system(command_bambu)

  expected_output_files <- file.path(outdir_test_data_bambu, c("counts_transcript.txt", "extended_annotations.gtf", "logfile.txt"))
  if (!base::all(file.exists(expected_output_files))) {
    error_msg <- "Error: Bambu failed to run. Some results files are missing. Please ensure all installation instructions are followed."
    return(error_msg)
  }

  # rename bambu output files
  renamed_gtf <- file.path(outdir_test_data_bambu, "bambu_transcript_annotations.gtf")
  file.rename(file.path(outdir_test_data_bambu,    "counts_transcript.txt"), file.path(outdir_test_data_bambu, "bambu_transcript_counts.txt"))
  file.rename(file.path(outdir_test_data_bambu, "extended_annotations.gtf"), renamed_gtf)

  # run gffcompare
  command_gff_compare <- paste0(conda_command, "gffcompare",
                                " -r ", shQuote(reference_transcriptome_file),
                                " ", shQuote(renamed_gtf))
  message(command_gff_compare)
  system(command_gff_compare)

  expected_output_file <- file.path(outdir_test_data_bambu, "gffcmp.bambu_transcript_annotations.gtf.tmap")
  if (!file.exists(expected_output_file)) {
    error_msg <- "Error: gffcompare failed to run. Some results files are missing. Please ensure all installation instructions are followed."
    return(error_msg)
  }

  # rename gffcompare output files
  file.rename(file.path(outdir_test_data_bambu, "gffcmp.bambu_transcript_annotations.gtf.tmap"), file.path(outdir_test_data_bambu, "gffcompare.tmap.txt"))
  file.remove(Sys.glob("gffcmp*"))

  ### database server ###
  outdir_test_data_db <- file.path(outdir_test_data, "database_output")
  if (!dir.exists(outdir_test_data_db)) {
    dir.create(outdir_test_data_db)
  }

  db_gtf_file <- file.path(outdir_test_data_bambu, "bambu_transcript_annotations.gtf")
  db_counts_file <- file.path(outdir_test_data_bambu, "bambu_transcript_counts.txt")

  ref_gtf <- reference_transcriptome_file
  counts_arg <- paste0(" -c ", shQuote(db_counts_file))
  vcf_arg <- paste0(" -v ", shQuote(vcf_file), " -G ", shQuote(reference_genome_file))

  command_generate_proteome <- paste0(conda_command, "Rscript bin/database_module/generate_proteome.R",
                                      " -g ", shQuote(db_gtf_file),
                                      " -r ", shQuote(ref_gtf),
                                      counts_arg,
                                      " -m ", shQuote(5),
                                      " -o ", shQuote("HUMAN"),
                                      " -l ", shQuote(30),
                                      " -u ", shQuote(FALSE),
                                      " -d ", shQuote(FALSE),
                                      vcf_arg,
                                      " -s ", shQuote(outdir_test_data_db))

  # run command
  message(command_generate_proteome)
  system(command_generate_proteome)

  expected_output_files <- file.path(outdir_test_data_db, c("ORFome_aa.txt", "proteome_database_transcripts.gtf", "Mutant_ORFome_aa.txt"))
  if (!base::all(file.exists(expected_output_files))) {
    error_msg <- "Error: Failed to generate proteome. Some results files are missing. Please ensure all installation instructions are followed."
    return(error_msg)
  }

  message("Generated ORFs")

  ref_proteome <- file.path("data", "openprot_uniprotDb_human.txt")
  orfome_file <- file.path(outdir_test_data_db, "ORFome_aa.txt")
  proteome_database_transcripts_gtf_file <- file.path(outdir_test_data_db, "proteome_database_transcripts.gtf")
  mutant_orfome_file <- file.path(outdir_test_data_db, "Mutant_ORFome_aa.txt")

  command_annotate_proteome <- paste0(conda_command, "python bin/database_module/annotate_proteome.py",
                                      " ", shQuote(ref_gtf),
                                      " ", shQuote(ref_proteome),
                                      " ", shQuote(orfome_file),
                                      " ", shQuote(proteome_database_transcripts_gtf_file),
                                      " ", shQuote(outdir_test_data_db),
                                      " ", shQuote("all"),
                                      " ", shQuote(30),
                                      " ", shQuote(mutant_orfome_file),
                                      " ", shQuote("HUMAN"),
                                      " ", shQuote(1),
                                      " ", shQuote("2000"))

  # annotate the proteome
  message(command_annotate_proteome)
  system(command_annotate_proteome)

  proteome_database_fasta_file <- file.path(outdir_test_data_db, "proteome_database.fasta")
  proteome_database_metadata_file <- file.path(outdir_test_data_db, "proteome_database_metadata.txt")
  expected_output_files <- c(proteome_database_fasta_file, proteome_database_metadata_file, proteome_database_transcripts_gtf_file)
  if (!base::all(file.exists(expected_output_files))) {
    error_msg <- "Error: Failed to annotate proteome. Some results files are missing. Please ensure all installation instructions are followed."
    return(error_msg)
  }

  message("Annotated proteome")

  ### integration server ###
  outdir_test_data_integ <- file.path(outdir_test_data, "integ_output")
  if (!dir.exists(outdir_test_data_integ)) {
    dir.create(outdir_test_data_integ)
  }

  command_map_peptides <- paste0(conda_command, "Rscript bin/integration_module/map_peptides_generate_outputs.R",
                                 " -p ", shQuote(reformatted_peptide_results_file),
                                 " -m ", shQuote(proteome_database_metadata_file),
                                 " -g ", shQuote(proteome_database_transcripts_gtf_file),
                                 " -s ", shQuote(outdir_test_data_integ))

  # map peptides and perform integration
  message(command_map_peptides)
  system(command_map_peptides)

  expected_output_files <- file.path(outdir_test_data_integ, c("peptides.bed12", "ORFs.bed12", "transcripts.bed12", "transcripts_and_ORFs_for_isovis.gtf", "peptide_info.tsv", "combined_annotations.gtf"))
  if (!base::all(file.exists(expected_output_files))) {
    error_msg <- "Error: Proteogenomics integration failed. Some results files are missing. Please ensure all installation instructions are followed."
    return(error_msg)
  }

  top_level_dir <- getwd()
  report_rmd_file <- file.path(top_level_dir, "bin", "integration_module", "integration_summary_report.Rmd")
  report_outdir <- file.path(top_level_dir, outdir_test_data_integ)

  command_create_report <- paste0(conda_command, "Rscript bin/integration_module/generate_integration_summary_report.R",
                                  " -i ", shQuote(report_rmd_file),
                                  " -o ", shQuote(report_outdir))

  # create the summary report
  message(command_create_report)
  system(command_create_report)

  # create a zip file with results
  files_to_zip_int <- c("summary_report.html", "peptide_info.tsv", "report_images/",
                        "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                        "peptides.bed12", "ORFs.bed12", "transcripts.bed12", "ncORF_stats.xlsx")

  # set the path to the ZIP file (in the session_id directory)
  zipfile_path_int <- file.path("..", "test_data_results.zip")

  # temp change the working dir to outdir_test_data_integ
  setwd(outdir_test_data_integ)

  if (base::all(file.exists(files_to_zip_int))) {
    if (file.exists(zipfile_path_int)) {
      file.remove(zipfile_path_int)
    }

    # zip files
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
  } else {
    error_msg <- "Error: Failed to generate proteogenomics integration summary report. Please ensure all installation instructions are followed."
  }

  # go back to starting dir
  setwd(top_level_dir)

  return(error_msg)
}

is_database_generation_input_valid <- function(input) {
  error_msg <- ""

  # The sequencing type selected must be either long-read or short-read
  sequencing_type <- input$sequencing_type
  if (is.null(sequencing_type) || !(sequencing_type %in% c("long-read", "short-read"))) {
    error_msg <- "Error: Please select the sequencing type (long-read or short-read)."
    return(error_msg)
  }

  # The input type selected must be FASTQs, BAMs or GTF
  input_type <- input$input_type
  if (is.null(input_type) || !(input_type %in% c("fastq_input", "bam_input", "gtf_input"))) {
    error_msg <- "Error: Please select the input type (FASTQs, BAMs or GTF)."
    return(error_msg)
  }

  # The organism selected must be human, roundworm, fruit fly, mouse, rat or zebrafish
  organism <- input$organism
  if (is.null(organism) || !(organism %in% c("HUMAN", "CAEEL", "DROME", "MOUSE", "RAT", "DANRE"))) {
    error_msg <- "Error: Please select the organism (human, roundworm, fruit fly, mouse, rat or zebrafish)."
    return(error_msg)
  }

  # Either canonical ORFs or all ORFs can be included in the proteome database
  database_type <- input$database_type
  if (is.null(database_type) || !(database_type %in% c("canonical", "all"))) {
    error_msg <- "Error: Please select the type of ORFs to include in the proteome database (canonical or all)."
    return(error_msg)
  }

  # Minimum ORF length: Integer >= 0
  min_orf_length <- floor(input$min_orf_length)
  if (!(is.finite(min_orf_length) && (min_orf_length >= 0))) {
    error_msg <- paste0("Error: The minimum ORF length must be a non-negative integer. Entered value: ", min_orf_length)
    return(error_msg)
  }

  # Minimum expression threshold: Integer >= 0
  minimum_tx_count <- floor(input$minimum_tx_count)
  if (!(is.finite(minimum_tx_count) && (minimum_tx_count >= 0))) {
    error_msg <- paste0("Error: The minimum expression threshold must be a non-negative integer. Entered value: ", minimum_tx_count)
    return(error_msg)
  }

  vcf_option <- input$vcf_option

  if (input_type == "fastq_input") {
    # Number of CPUs: Positive integer
    user_threads <- floor(input$user_threads)
    if (!(is.finite(user_threads) && (user_threads >= 1))) {
      error_msg <- paste0("Error: The number of CPUs must be a positive integer. Entered value: ", user_threads)
      return(error_msg)
    }

    # Check if the user has uploaded a reference genome FASTA file yet
    user_reference_genome <- input$user_reference_genome
    if (is.null(user_reference_genome)) {
      error_msg <- "Error: Please upload a reference genome FASTA file for the proteome database generation."
      return(error_msg)
    }

    # The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped)
    filename <- user_reference_genome$name
    if (!any(endsWith(filename, c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")))) {
      error_msg <- "Error: The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped."
      return(error_msg)
    }

    # The reference genome FASTA file must not be empty
    file_path <- user_reference_genome$datapath
    if (file.size(file_path) == 0) {
      error_msg <- "Error: The reference genome FASTA file must not be empty."
      return(error_msg)
    }

    if (sequencing_type == "short-read") {
      # Check if the user has uploaded a reference transcriptome FASTA file yet
      transcriptome_file <- input$transcriptome_file
      if (is.null(transcriptome_file)) {
        error_msg <- "Error: Please upload a reference transcriptome FASTA file for the proteome database generation."
        return(error_msg)
      }

      # The reference transcriptome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped)
      filename <- transcriptome_file$name
      if (!any(endsWith(filename, c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")))) {
        error_msg <- "Error: The reference transcriptome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped."
        return(error_msg)
      }

      # The reference transcriptome FASTA file must not be empty
      file_path <- transcriptome_file$datapath
      if (file.size(file_path) == 0) {
        error_msg <- "Error: The reference transcriptome FASTA file must not be empty."
        return(error_msg)
      }
    }

    # Check if the user has uploaded FASTQ files yet
    user_fastq_files <- input$user_fastq_files
    if (is.null(user_fastq_files)) {
      error_msg <- "Error: At least one FASTQ file must be uploaded."
      return(error_msg)
    }

    # Every FASTQ file must have a file extension of '.fastq', '.fq', '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped), and must not be empty
    num_files <- nrow(user_fastq_files)
    filenames <- user_fastq_files$name
    for (i in 1:num_files) {
      filename <- filenames[i]
      if (!any(endsWith(filename, c(".fastq", ".fq", ".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn", ".fastq.gz", ".fq.gz", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")))) {
        error_msg <- "Error: Please ensure that all uploaded FASTQ files have a file extension of '.fastq', '.fq', '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped."
        return(error_msg)
      }

      file_path <- user_fastq_files$datapath[i]
      if (file.size(file_path) == 0) {
        error_msg <- "Error: Please ensure that none of the uploaded FASTQ files are empty."
        return(error_msg)
      }
    }
  } else if (input_type == "bam_input") {
    # Number of CPUs: Positive integer
    user_threads <- floor(input$user_threads)
    if (!(is.finite(user_threads) && (user_threads >= 1))) {
      error_msg <- paste0("Error: The number of CPUs must be a positive integer. Entered value: ", user_threads)
      return(error_msg)
    }

    # Check if the user has uploaded a reference genome FASTA file yet
    user_reference_genome <- input$user_reference_genome
    if (is.null(user_reference_genome)) {
      error_msg <- "Error: Please upload a reference genome FASTA file for the proteome database generation."
      return(error_msg)
    }

    # The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped)
    filename <- user_reference_genome$name
    if (!any(endsWith(filename, c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")))) {
      error_msg <- "Error: The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped."
      return(error_msg)
    }

    # The reference genome FASTA file must not be empty
    file_path <- user_reference_genome$datapath
    if (file.size(file_path) == 0) {
      error_msg <- "Error: The reference genome FASTA file must not be empty."
      return(error_msg)
    }

    # Check if the user has uploaded BAM files yet
    user_bam_files <- input$user_bam_files
    if (is.null(user_bam_files)) {
      error_msg <- "Error: At least one BAM file must be uploaded."
      return(error_msg)
    }

    # Every BAM file must have a file extension of '.bam' and must not be empty
    num_files <- nrow(user_bam_files)
    filenames <- user_bam_files$name
    for (i in 1:num_files) {
      filename <- filenames[i]
      if (!endsWith(filename, ".bam")) {
        error_msg <- "Error: Please ensure that all uploaded BAM files have a file extension of '.bam'."
        return(error_msg)
      }

      file_path <- user_bam_files$datapath[i]
      if (file.size(file_path) == 0) {
        error_msg <- "Error: Please ensure that none of the uploaded BAM files are empty."
        return(error_msg)
      }
    }
  } else {
    if (sequencing_type == "long-read") {
      # Check if the user has uploaded a user-generated transcript annotation GTF file yet
      user_gtf_file <- input$user_gtf_file
      if (is.null(user_gtf_file)) {
        error_msg <- "Error: Please upload a user-generated transcript annotation GTF file for the proteome database generation."
        return(error_msg)
      }

      # The user-generated transcript annotation GTF file must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'
      filename <- user_gtf_file$name
      if (!any(endsWith(filename, c(".gtf", ".gff", ".gff2", ".gff3")))) {
        error_msg <- "Error: The user-generated transcript annotation GTF file must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'."
        return(error_msg)
      }

      # The user-generated transcript annotation GTF file must not be empty
      file_path <- user_gtf_file$datapath
      if (file.size(file_path) == 0) {
        error_msg <- "Error: The user-generated transcript annotation GTF file must not be empty."
        return(error_msg)
      }

      # If the user has uploaded a user-generated transcript counts file...
      user_tx_count_file <- input$user_tx_count_file
      if (!is.null(user_tx_count_file)) {
        # The user-generated transcript counts file must have a file extension of '.txt', '.csv' or '.tsv'
        filename <- user_tx_count_file$name
        if (!any(endsWith(filename, c(".txt", ".csv", ".tsv")))) {
          error_msg <- "Error: The user-generated transcript counts file must have a file extension of '.txt', '.csv' or '.tsv'."
          return(error_msg)
        }

        # The user-generated transcript counts file must not be empty
        file_path <- user_tx_count_file$datapath
        if (file.size(file_path) == 0) {
          error_msg <- "Error: The user-generated transcript counts file must not be empty."
          return(error_msg)
        }
      }
    } else {
      # Check if the user has uploaded a user-generated transcript counts file yet
      user_tx_count_file <- input$user_tx_count_file
      if (is.null(user_tx_count_file)) {
        error_msg <- "Error: Please upload a user-generated transcript counts file for the proteome database generation."
        return(error_msg)
      }

      # The user-generated transcript counts file must have a file extension of '.txt', '.csv' or '.tsv'
      filename <- user_tx_count_file$name
      if (!any(endsWith(filename, c(".txt", ".csv", ".tsv")))) {
        error_msg <- "Error: The user-generated transcript counts file must have a file extension of '.txt', '.csv' or '.tsv'."
        return(error_msg)
      }

      # The user-generated transcript counts file must not be empty
      file_path <- user_tx_count_file$datapath
      if (file.size(file_path) == 0) {
        error_msg <- "Error: The user-generated transcript counts file must not be empty."
        return(error_msg)
      }
    }

    if (vcf_option) {
      # Check if the user has uploaded a reference genome FASTA file yet
      user_reference_genome <- input$user_reference_genome
      if (is.null(user_reference_genome)) {
        error_msg <- "Error: Please upload a reference genome FASTA file for the proteome database generation."
        return(error_msg)
      }

      # The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped)
      filename <- user_reference_genome$name
      if (!any(endsWith(filename, c(".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna.gz", ".ffn.gz", ".faa.gz", ".mpfa.gz", ".frn.gz")))) {
        error_msg <- "Error: The reference genome FASTA file must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped."
        return(error_msg)
      }

      # The reference genome FASTA file must not be empty
      file_path <- user_reference_genome$datapath
      if (file.size(file_path) == 0) {
        error_msg <- "Error: The reference genome FASTA file must not be empty."
        return(error_msg)
      }
    }
  }

  # Check if the user has uploaded a reference transcriptome GTF file yet
  reference_gtf_file <- input$reference_gtf_file
  if (is.null(reference_gtf_file)) {
    error_msg <- "Error: Please upload a reference transcriptome GTF file for the proteome database generation."
    return(error_msg)
  }

  # The reference transcriptome GTF file must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'
  filename <- reference_gtf_file$name
  if (!any(endsWith(filename, c(".gtf", ".gff", ".gff2", ".gff3")))) {
    error_msg <- "Error: The reference transcriptome GTF file must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'."
    return(error_msg)
  }

  # The reference transcriptome GTF file must not be empty
  file_path <- reference_gtf_file$datapath
  if (file.size(file_path) == 0) {
    error_msg <- "Error: The reference transcriptome GTF file must not be empty."
    return(error_msg)
  }

  # If the user has selected to upload a VCF file...
  if (vcf_option) {
    # Check if the user has uploaded a VCF file yet
    user_vcf_file <- input$user_vcf_file
    if (is.null(user_vcf_file)) {
      error_msg <- "Error: Please upload a VCF file for the proteome database generation."
      return(error_msg)
    }

    # The VCF file must have a file extension of '.vcf'
    filename <- user_vcf_file$name
    if (!endsWith(filename, ".vcf")) {
      error_msg <- "Error: The VCF file must have a file extension of '.vcf'."
      return(error_msg)
    }

    # The VCF file must not be empty
    file_path <- user_vcf_file$datapath
    if (file.size(file_path) == 0) {
      error_msg <- "Error: The VCF file must not be empty."
      return(error_msg)
    }
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

fastq_server <- function(input, session) {
  session_id <- session$token                           # session ID
  outdir_bam <- file.path(session_id, "mapping_output") # output directory

  if (!dir.exists(outdir_bam)) {    # create the output directory if it does not exist yet
    dir.create(outdir_bam)
  } else {                          # if it already exists, remove all Salmon quant.sf files inside it
    quant_sf_files <- Sys.glob(file.path(outdir_bam, "*", "quant.sf"))
    for (quant_sf_file in quant_sf_files) {
      if (file.exists(quant_sf_file)) {
        file.remove(quant_sf_file)
      }
    }
  }

  user_fastq_files_df <- input$user_fastq_files %>%
    dplyr::mutate(file_prefix = stringr::str_replace_all(name, c("\\.fastq$" = "", "\\.fq$" = "", "\\.fasta$" = "", "\\.fas$" = "", "\\.fa$" = "", "\\.fna$" = "", "\\.ffn$" = "", "\\.faa$" = "", "\\.mpfa$" = "", "\\.frn$" = "",
                                                                 "\\.fastq\\.gz$" = "", "\\.fq\\.gz$" = "", "\\.fasta\\.gz$" = "", "\\.fas\\.gz$" = "", "\\.fa\\.gz$" = "", "\\.fna\\.gz$" = "", "\\.ffn\\.gz$" = "", "\\.faa\\.gz$" = "", "\\.mpfa\\.gz$" = "", "\\.frn\\.gz$" = "")))

  print(user_fastq_files_df)

  # map reads
  user_threads <- floor(input$user_threads)
  user_reference_genome <- input$user_reference_genome
  reference_gtf_file <- input$reference_gtf_file$datapath
  transcriptome_file <- input$transcriptome_file$datapath
  if (input$sequencing_type == "long-read") {
    genome_index_file <- file.path(outdir_bam, paste0(stringr::str_replace_all(user_reference_genome$name, c("\\.fasta$" = "", "\\.fas$" = "", "\\.fa$" = "", "\\.fna$" = "", "\\.ffn$" = "", "\\.faa$" = "", "\\.mpfa$" = "", "\\.frn$" = "",
                                                                                                             "\\.fasta\\.gz$" = "", "\\.fas\\.gz$" = "", "\\.fa\\.gz$" = "", "\\.fna\\.gz$" = "", "\\.ffn\\.gz$" = "", "\\.faa\\.gz$" = "", "\\.mpfa\\.gz$" = "", "\\.frn\\.gz$" = "")),
                                                      ".mmi"))
    command_minimap2_index <- paste0(conda_command, "minimap2",
                                     " -ax ", shQuote("splice:hq"),
                                     " -d ", shQuote(genome_index_file),
                                     " ", shQuote(user_reference_genome$datapath))

    # index the genome
    message(command_minimap2_index)
    system(command_minimap2_index)

    # for each FASTQ file, go through the following 3 steps...
    for (i in 1:nrow(user_fastq_files_df)) {
      fastq_file <- user_fastq_files_df$datapath[i]
      file_prefix <- user_fastq_files_df$file_prefix[i]

      # 1. create a SAM file of reads mapped to the indexed genome
      mapped_reads_file <- file.path(outdir_bam, paste0(file_prefix, ".sam"))
      command_minimap2 <- paste0(conda_command, "minimap2 --sam-hit-only --secondary=no",
                                 " -t ", shQuote(user_threads),
                                 " -ax ", shQuote("splice:hq"),
                                 " -o ", shQuote(mapped_reads_file),
                                 " ", shQuote(genome_index_file),
                                 " ", shQuote(fastq_file))
      message(command_minimap2)
      system(command_minimap2)

      # 2. turn the SAM file into a BAM file with headers
      unsorted_bam_file <- file.path(outdir_bam, paste0(file_prefix, "_unsorted.bam"))
      command_samtools_view <- paste0(conda_command, "samtools view -bh",
                                      " -o ", shQuote(unsorted_bam_file),
                                      " ", shQuote(mapped_reads_file))
      message(command_samtools_view)
      system(command_samtools_view)
      if (file.exists(mapped_reads_file)) {
        file.remove(mapped_reads_file)
      }

      # 3. sort the reads in the BAM file by ascending genomic coordinates
      sorted_bam_file <- file.path(outdir_bam, paste0(file_prefix, ".bam"))
      command_samtools_sort <- paste0(conda_command, "samtools sort",
                                      " -@ ", shQuote(user_threads),
                                      " -o ", shQuote(sorted_bam_file),
                                      " ", shQuote(unsorted_bam_file))
      message(command_samtools_sort)
      system(command_samtools_sort)
      if (file.exists(unsorted_bam_file)) {
        file.remove(unsorted_bam_file)
      }
    }
  } else {
    is_reference_genome_gzipped <- grepl("\\.gz$", user_reference_genome$datapath)
    is_reference_transcriptome_gzipped <- grepl("\\.gz$", transcriptome_file)
    decoy_file <- file.path(outdir_bam, "decoys.txt")
    salmon_index_file <- file.path(outdir_bam, "salmon_index")

    if (is_reference_genome_gzipped) {  # if the genome is gzipped, decompress it before generating decoys
      temp_decompressed_genome_file <- file.path(outdir_bam, "temp_genome.fa")
      command_decompress <- paste0("gzip -c -d ", shQuote(user_reference_genome$datapath), " > ", shQuote(temp_decompressed_genome_file))
      message(command_decompress)
      system(command_decompress)

      command_generate_decoy <- paste0("grep '^>' ",  shQuote(temp_decompressed_genome_file), " | cut -d ' ' -f 1 > ", shQuote(decoy_file))
    } else {
      command_generate_decoy <- paste0("grep '^>' ", shQuote(user_reference_genome$datapath), " | cut -d ' ' -f 1 > ", shQuote(decoy_file))
    }

    command_sed <- paste0("sed -i -e 's/>//g' ", shQuote(decoy_file))

    genome_transcriptome_combined_file <- file.path(outdir_bam, "gentrome.fa")
    if (is_reference_genome_gzipped & is_reference_transcriptome_gzipped) {             # if both the genome and transcriptome are gzipped, combine them directly
      genome_transcriptome_combined_file <- file.path(outdir_bam, "gentrome.fa.gz")

      command_ref_file <- paste0("cat ",                   shQuote(transcriptome_file), " ", shQuote(user_reference_genome$datapath), " > ", shQuote(genome_transcriptome_combined_file))
    } else if (!is_reference_genome_gzipped & !is_reference_transcriptome_gzipped) {    # if both the genome and transcriptome are not gzipped, combine them directly
      command_ref_file <- paste0("cat ",                   shQuote(transcriptome_file), " ", shQuote(user_reference_genome$datapath), " > ", shQuote(genome_transcriptome_combined_file))
    } else if (is_reference_genome_gzipped) {                                           # if only the genome is gzipped, decompress it before combining it with the transcriptome
      command_ref_file <- paste0("cat ",                   shQuote(transcriptome_file), " ",  shQuote(temp_decompressed_genome_file), " > ", shQuote(genome_transcriptome_combined_file))
    } else {                                                                            # if only the transcriptome is gzipped, decompress it before combining it with the genome
      temp_decompressed_transcriptome_file <- file.path(outdir_bam, "temp_transcriptome.fa")
      command_decompress <- paste0("gzip -c -d ", shQuote(transcriptome_file), " > ", shQuote(temp_decompressed_transcriptome_file))
      message(command_decompress)
      system(command_decompress)

      command_ref_file <- paste0("cat ", shQuote(temp_decompressed_transcriptome_file), " ", shQuote(user_reference_genome$datapath), " > ", shQuote(genome_transcriptome_combined_file))
    }

    command_index <- paste0(conda_command, "salmon index --gencode",
                            " -t ", shQuote(genome_transcriptome_combined_file),
                            " -d ", shQuote(decoy_file),
                            " -p ", shQuote(user_threads),
                            " -i ", shQuote(salmon_index_file))

    # generate the salmon index
    commands <- c(command_generate_decoy, command_sed, command_ref_file, command_index)
    for (i in 1:length(commands)) {
      message(commands[i])
      system(commands[i])

      # remove the temporary decompressed genome or transcriptome file after creating the combined genome and transcriptome file
      if (i == 3) {
        if (is_reference_genome_gzipped & file.exists(temp_decompressed_genome_file)) {
          file.remove(temp_decompressed_genome_file)
        } else if (!is_reference_genome_gzipped & is_reference_transcriptome_gzipped & file.exists(temp_decompressed_transcriptome_file)) {
          file.remove(temp_decompressed_transcriptome_file)
        }
      }
    }

    # for each FASTQ file, determine whether it contains paired-end reads or single-end reads
    # FIXME: handle FASTQ files with no sample name (e.g. ".fastq", "_R1.fastq", "_R2.fastq")
    # FIXME: handle samples that contain both paired-end reads and single-end reads
    paired_end <- list()
    single_end <- list()

    for (i in 1:nrow(user_fastq_files_df)) {
      file_name <- user_fastq_files_df$name[i]
      if (grepl("_R1", file_name)) {
        base_name <- stringr::str_replace_all(file_name, c("_R1\\.fastq$" = "", "_R1\\.fq$" = "", "_R1\\.fasta$" = "", "_R1\\.fas$" = "", "_R1\\.fa$" = "", "_R1\\.fna$" = "", "_R1\\.ffn$" = "", "_R1\\.faa$" = "", "_R1\\.mpfa$" = "", "_R1\\.frn$" = "",
                                                           "_R1\\.fastq\\.gz$" = "", "_R1\\.fq\\.gz$" = "", "_R1\\.fasta\\.gz$" = "", "_R1\\.fas\\.gz$" = "", "_R1\\.fa\\.gz$" = "", "_R1\\.fna\\.gz$" = "", "_R1\\.ffn\\.gz$" = "", "_R1\\.faa\\.gz$" = "", "_R1\\.mpfa\\.gz$" = "", "_R1\\.frn\\.gz$" = ""))
        paired_end[[base_name]]$R1 <- user_fastq_files_df$datapath[i]
      } else if (grepl("_R2", file_name)) {
        base_name <- stringr::str_replace_all(file_name, c("_R2\\.fastq$" = "", "_R2\\.fq$" = "", "_R2\\.fasta$" = "", "_R2\\.fas$" = "", "_R2\\.fa$" = "", "_R2\\.fna$" = "", "_R2\\.ffn$" = "", "_R2\\.faa$" = "", "_R2\\.mpfa$" = "", "_R2\\.frn$" = "",
                                                           "_R2\\.fastq\\.gz$" = "", "_R2\\.fq\\.gz$" = "", "_R2\\.fasta\\.gz$" = "", "_R2\\.fas\\.gz$" = "", "_R2\\.fa\\.gz$" = "", "_R2\\.fna\\.gz$" = "", "_R2\\.ffn\\.gz$" = "", "_R2\\.faa\\.gz$" = "", "_R2\\.mpfa\\.gz$" = "", "_R2\\.frn\\.gz$" = ""))
        paired_end[[base_name]]$R2 <- user_fastq_files_df$datapath[i]
      } else {
        base_name <- user_fastq_files_df$file_prefix[i]
        single_end[[base_name]] <- user_fastq_files_df$datapath[i]
      }
    }

    # quantify paired-end reads
    for (base_name in names(paired_end)) {
      if (base_name == "") {
        next
      }

      R1_path <- paired_end[[base_name]]$R1
      R2_path <- paired_end[[base_name]]$R2

      if (!is.null(R1_path) & !is.null(R2_path)) {
        quant_output_folder <- file.path(outdir_bam, base_name)
        command_salmon <- paste0(conda_command, "salmon quant --validateMappings",
                                 " -i ", shQuote(salmon_index_file),
                                 " -p ", shQuote(user_threads),
                                 " -l ", shQuote("A"),
                                 " -1 ", shQuote(R1_path),
                                 " -2 ", shQuote(R2_path),
                                 " -o ", shQuote(quant_output_folder))

        message(paste0("Base name: ", base_name, " --- ", "R1 path: ", R1_path, ", R2 path: ", R2_path))
        message(command_salmon)
        system(command_salmon)
      }
    }

    # quantify single-end reads
    for (base_name in names(single_end)) {
      if (base_name == "") { 
        next
      }

      quant_output_folder <- file.path(outdir_bam, base_name)
      command_salmon <- paste0(conda_command, "salmon quant --validateMappings",
                               " -i ", shQuote(salmon_index_file),
                               " -p ", shQuote(user_threads),
                               " -l ", shQuote("A"),
                               " -r ", shQuote(single_end[[base_name]]),
                               " -o ", shQuote(quant_output_folder))

      message(single_end[[base_name]])
      message(command_salmon)
      system(command_salmon)
    }

    transcript_counts_file <- file.path(outdir_bam, "bambu_transcript_counts.txt")
    command_create_count_matrix <- paste0(conda_command, "Rscript bin/database_module/matrix_compilation_salmon.R",
                                          " -s ", shQuote(outdir_bam),
                                          " -g ", shQuote(reference_gtf_file),
                                          " -o ", shQuote(transcript_counts_file))

    # create the count matrix
    message(command_create_count_matrix)
    system(command_create_count_matrix)
  }
}

bam_server <- function(input, session) {
  session_id <- session$token                           # session ID
  outdir_bam <- file.path(session_id, "mapping_output") # output directory

  if (!dir.exists(outdir_bam)) {    # create the output directory if it does not exist yet
    dir.create(outdir_bam)
  } else {                          # if it already exists, remove all Salmon quant.sf files inside it
    quant_sf_files <- Sys.glob(file.path(outdir_bam, "*", "quant.sf"))
    for (quant_sf_file in quant_sf_files) {
      if (file.exists(quant_sf_file)) {
        file.remove(quant_sf_file)
      }
    }
  }

  # generate the reference transcript FASTA file
  user_threads <- floor(input$user_threads)
  user_reference_genome <- input$user_reference_genome$datapath
  reference_gtf_file <- input$reference_gtf_file$datapath
  transcript_fasta_file <- file.path(outdir_bam, "transcript.fa")

  is_reference_genome_gzipped <- grepl("\\.gz$", user_reference_genome)
  if (is_reference_genome_gzipped) {    # if the reference genome file is gzipped, decompress it first
    decompressed_genome_fasta_file <- file.path(outdir_bam, "temp_genome.fa")
    command_decompress <- paste0("gzip -c -d ", shQuote(user_reference_genome), " > ", shQuote(decompressed_genome_fasta_file))
    message(command_decompress)
    system(command_decompress)

    command_gffread <- paste0(conda_command, "gffread",
                              " -w ", shQuote(transcript_fasta_file),
                              " -g ", shQuote(decompressed_genome_fasta_file),
                              " ", shQuote(reference_gtf_file))
  } else {
    command_gffread <- paste0(conda_command, "gffread",
                              " -w ", shQuote(transcript_fasta_file),
                              " -g ", shQuote(user_reference_genome),
                              " ", shQuote(reference_gtf_file))
  }

  message(command_gffread)
  system(command_gffread)

  # remove the temporary decompressed reference genome file
  if (is_reference_genome_gzipped & file.exists(decompressed_genome_fasta_file)) {
    file.remove(decompressed_genome_fasta_file)
  }

  # create df of bam file names
  user_bam_files_df <- input$user_bam_files %>%
    dplyr::mutate(file_prefix = sub("\\.bam$", "", name))

  # for each bam file, quantify reads
  for (i in 1:nrow(user_bam_files_df)) {
    bam_file <- user_bam_files_df$datapath[i]
    file_prefix <- user_bam_files_df$file_prefix[i]
    quant_output_folder <- file.path(outdir_bam, file_prefix)

    command_salmon <- paste0(conda_command, "salmon quant",
                             " -t ", shQuote(transcript_fasta_file),
                             " -p ", shQuote(user_threads),
                             " -l ", shQuote("A"),
                             " -a ", shQuote(bam_file),
                             " -o ", shQuote(quant_output_folder))

    message(command_salmon)
    system(command_salmon)
  }

  transcript_counts_file <- file.path(outdir_bam, "bambu_transcript_counts.txt")
  command_create_count_matrix <- paste0(conda_command, "Rscript bin/database_module/matrix_compilation_salmon.R",
                                        " -s ", shQuote(outdir_bam),
                                        " -g ", shQuote(reference_gtf_file),
                                        " -o ", shQuote(transcript_counts_file))

  # create the count matrix
  message(command_create_count_matrix)
  system(command_create_count_matrix)
}

bambu_server <- function(input, session) {
  session_id <- session$token                           # session ID
  outdir_bam <- file.path(session_id, "mapping_output")
  outdir_bambu <- file.path(session_id, "bambu_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_bambu)) {
    dir.create(outdir_bambu)
  }

  user_threads <- floor(input$user_threads)
  user_bam_files <- input$user_bam_files
  reference_gtf_file <- input$reference_gtf_file$datapath
  organism <- input$organism
  logfile_path <- file.path(outdir_bambu, "logfile.txt")

  if (input$input_type == "bam_input") {    # if the user uploaded BAM files
    bamdir <- dirname(user_bam_files$datapath)

    # Specify new filenames
    new_names <- file.path(bamdir[1], user_bam_files$name)

    # Rename the files
    file.rename(user_bam_files$datapath, new_names)

    command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R",
                            " -b ", shQuote(bamdir[1]),
                            " -g ", shQuote(reference_gtf_file),
                            " -o ", shQuote(outdir_bambu),
                            " -t ", shQuote(user_threads),
                            " -s ", shQuote(organism),
                            " > ", shQuote(logfile_path),
                            " 2>&1")
  } else {                                  # if the user uploaded FASTQ files
    command_bambu <- paste0(conda_command, "Rscript bin/database_module/run_bambu.R",
                            " -b ", shQuote(outdir_bam),
                            " -g ", shQuote(reference_gtf_file),
                            " -o ", shQuote(outdir_bambu),
                            " -t ", shQuote(user_threads),
                            " -s ", shQuote(organism),
                            " > ", shQuote(logfile_path),
                            " 2>&1")
  }

  # run bambu
  message(command_bambu)
  system(command_bambu)

  # rename bambu output files
  renamed_gtf <- file.path(outdir_bambu, "bambu_transcript_annotations.gtf")
  file.rename(file.path(outdir_bambu,    "counts_transcript.txt"), file.path(outdir_bambu, "bambu_transcript_counts.txt"))
  file.rename(file.path(outdir_bambu, "extended_annotations.gtf"), renamed_gtf)

  # run gffcompare
  command_gff_compare <- paste0(conda_command, "gffcompare",
                                " -r ", shQuote(reference_gtf_file),
                                " ", shQuote(renamed_gtf))
  message(command_gff_compare)
  system(command_gff_compare)

  # rename gffcompare output files
  file.rename(file.path(outdir_bambu, "gffcmp.bambu_transcript_annotations.gtf.tmap"), file.path(outdir_bambu, "gffcompare.tmap.txt"))
  file.remove(Sys.glob("gffcmp*"))
}

database_server <- function(input, session) {
  session_id <- session$token                           # session ID
  outdir_db <- file.path(session_id, "database_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_db)) {
    dir.create(outdir_db)
  }

  if (input$input_type == "gtf_input" && input$sequencing_type == "long-read") {    # if the user supplied GTF files and long-read RNA-seq data
    db_gtf_file <- input$user_gtf_file$datapath
    db_counts_file <- input$user_tx_count_file$datapath
  } else if (input$sequencing_type == "long-read") {                                # if the user supplied FASTQ / BAM files and long-read RNA-seq data
    db_gtf_file <- file.path(session_id, "bambu_output", "bambu_transcript_annotations.gtf")
    db_counts_file <- file.path(session_id, "bambu_output", "bambu_transcript_counts.txt")
  } else {                                                                          # if the user supplied short-read RNA-seq data
    db_gtf_file <- input$reference_gtf_file$datapath
    db_counts_file <- file.path(session_id, "mapping_output", "bambu_transcript_counts.txt")
  }

  ref_gtf <- input$reference_gtf_file$datapath

  # conditionally add the transcript counts file argument
  counts_arg <- if (!is.null(db_counts_file) && nzchar(db_counts_file)) {
    paste0(" -c ", shQuote(db_counts_file))
  } else {
    ""
  }

  # conditionally add the VCF file argument (the reference genome will only be used if a VCF file was provided)
  vcf_arg <- if (input$vcf_option && !is.null(input$user_vcf_file$datapath) && !is.null(input$user_reference_genome$datapath)) {
    vcf_file <- input$user_vcf_file$datapath

    is_reference_genome_gzipped <- grepl("\\.gz$", input$user_reference_genome$datapath)
    if (is_reference_genome_gzipped) {  # if the genome is gzipped, decompress it
      temp_decompressed_genome_file <- file.path(outdir_db, "temp_genome.fa")
      command_decompress <- paste0("gzip -c -d ", shQuote(input$user_reference_genome$datapath), " > ", shQuote(temp_decompressed_genome_file))
      message(command_decompress)
      system(command_decompress)

      paste0(" -v ", shQuote(vcf_file), " -G ",        shQuote(temp_decompressed_genome_file))
    } else {
      paste0(" -v ", shQuote(vcf_file), " -G ", shQuote(input$user_reference_genome$datapath))
    }
  } else {
    vcf_file <- NULL
    ""
  }

  # construct the command
  minimum_tx_count <- floor(input$minimum_tx_count)
  min_orf_length <- floor(input$min_orf_length)
  organism <- input$organism
  user_threads <- floor(input$user_threads)

  command_generate_proteome <- paste0(conda_command, "Rscript bin/database_module/generate_proteome.R",
                                      " -g ", shQuote(db_gtf_file),
                                      " -r ", shQuote(ref_gtf),
                                      counts_arg,
                                      " -m ", shQuote(minimum_tx_count),
                                      " -o ", shQuote(organism),
                                      " -l ", shQuote(min_orf_length),
                                      " -u ", shQuote(input$user_find_utr_5_orfs),
                                      " -d ", shQuote(input$user_find_utr_3_orfs),
                                      vcf_arg,
                                      " -s ", shQuote(outdir_db))

  # run command
  message(command_generate_proteome)
  system(command_generate_proteome)
  message("Generated ORFs")

  # set reference protein database per organism
  if (organism == "HUMAN") {
    ref_proteome <- file.path("data", "openprot_uniprotDb_human.txt")
  } else if (organism == "CAEEL") {
    ref_proteome <- file.path("data", "openprot_uniprotDb_c_elegans.txt")
  } else if (organism == "DROME") {
    ref_proteome <- file.path("data", "openprot_uniprotDb_drosophila.txt")
  } else if (organism == "MOUSE") {
    ref_proteome <- file.path("data", "openprot_uniprotDb_mouse.txt")
  } else if (organism == "RAT") {
    ref_proteome <- file.path("data", "openprot_uniprotDb_rat.txt")
  } else {
    ref_proteome <- file.path("data", "openprot_uniprotDb_zebrafish.txt")
  }
  #} else if (organism == "DANRE") {
  #  ref_proteome <- file.path("data", "openprot_uniprotDb_zebrafish.txt")
  #} else if (organism == "PANTR") {
  #  ref_proteome <- file.path("data", "openprot_uniprotDb_chimp.txt")
  #} else if (organism == "BOVIN") {
  #  ref_proteome <- file.path("data", "openprot_uniprotDb_cow.txt")
  #} else if (organism == "XENTR") {
  #  ref_proteome <- file.path("data", "openprot_uniprotDb_clawed_frog.txt")
  #} else {
  #  ref_proteome <- file.path("data", "openprot_uniprotDb_yeast.txt")
  #}

  orfome_file <- file.path(outdir_db, "ORFome_aa.txt")
  proteome_database_transcripts_gtf_file <- file.path(outdir_db, "proteome_database_transcripts.gtf")
  mutant_orfome_file <- file.path(outdir_db, "Mutant_ORFome_aa.txt")
  if (!is.null(vcf_file)) { # if there is a VCF file uploaded
    command_annotate_proteome <- paste0(conda_command, "python bin/database_module/annotate_proteome.py",
                                        " ", shQuote(ref_gtf),
                                        " ", shQuote(ref_proteome),
                                        " ", shQuote(orfome_file),
                                        " ", shQuote(proteome_database_transcripts_gtf_file),
                                        " ", shQuote(outdir_db),
                                        " ", shQuote(input$database_type),
                                        " ", shQuote(min_orf_length),
                                        " ", shQuote(mutant_orfome_file),
                                        " ", shQuote(organism),
                                        " ", shQuote(user_threads),
                                        " ", shQuote("2000"))
  } else { # if no VCF file uploaded
    command_annotate_proteome <- paste0(conda_command, "python bin/database_module/annotate_proteome.py",
                                        " ", shQuote(ref_gtf),
                                        " ", shQuote(ref_proteome),
                                        " ", shQuote(orfome_file),
                                        " ", shQuote(proteome_database_transcripts_gtf_file),
                                        " ", shQuote(outdir_db),
                                        " ", shQuote(input$database_type),
                                        " ", shQuote(min_orf_length),
                                        " ", shQuote("None"),
                                        " ", shQuote(organism),
                                        " ", shQuote(user_threads),
                                        " ", shQuote("2000"))
  }

  # annotate the proteome
  message(command_annotate_proteome)
  system(command_annotate_proteome)
  message("Annotated proteome")

  # remove the temporary decompressed genome file if a VCF file and a gzipped reference genome file were provided
  if (nchar(vcf_arg) != 0) {
    is_reference_genome_gzipped <- grepl("\\.gz$", input$user_reference_genome$datapath)
    if (is_reference_genome_gzipped) {
      temp_decompressed_genome_file <- file.path(outdir_db, "temp_genome.fa")
      if (file.exists(temp_decompressed_genome_file)) {
        file.remove(temp_decompressed_genome_file)
      }
    }
  }

  # get top level directory
  top_level_dir <- getwd()

  # zip all results files depending on input types
  proteome_database_fasta_file <- file.path(outdir_db, "proteome_database.fasta")
  orf_temp_file <- file.path(outdir_db, "orf_temp.txt")
  if (file.exists(proteome_database_fasta_file) & file.exists(proteome_database_transcripts_gtf_file) & !file.exists(orf_temp_file)) {
    files_to_zip_db <- c("proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")

    if (input$input_type == "fastq_input" & input$sequencing_type == "long-read") {
      bam_files <- list.files(path = file.path("..", "mapping_output"), "\\.bam$", full.names = TRUE)
      files_to_zip_db <- c(files_to_zip_db, bam_files, "../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/bambu_transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "../bambu_output/logfile.txt")
    } else if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      files_to_zip_db <- c(files_to_zip_db, "../bambu_output/bambu_transcript_annotations.gtf", "../bambu_output/bambu_transcript_counts.txt", "../bambu_output/novel_transcript_classes.csv", "../bambu_output/gffcompare.tmap.txt", "../bambu_output/logfile.txt")
    } else if (input$sequencing_type == "short-read") {
      files_to_zip_db <- c(files_to_zip_db, "../mapping_output/bambu_transcript_counts.txt")
    }

    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_db <- file.path("..", "database_results.zip")

    # temp change the working dir to outdir_db
    setwd(outdir_db)

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

  # The proteome database must not be empty
  file_path <- fragpipe_prot_db_file$datapath
  if (file.size(file_path) == 0) {
    error_msg <- "Error: The proteome database file must not be empty."
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

    # Mass spectrometry data files must not be empty
    file_path <- mass_spec_file$datapath
    if (file.size(file_path) == 0) {
      error_msg <- "Error: Please ensure that none of the uploaded mass spectrometry data files are empty."
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

fragpipe_server <- function(input, session_id, mass_spec_file_num) {
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
  output_dir <- "fragpipe_output"   # output directory

  # Create the output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Go into the output directory so that the FragPipe log file will be output there
  output_dir <- normalizePath(output_dir)
  setwd(output_dir)

  # Write the mass spectrometry file list
  mass_spec_info_file <- "mass_spec_info_list.txt"
  writeLines(mass_spec_info_file_contents, mass_spec_info_file)

  command_run_fragpipe <- paste0(conda_command, "python ../../bin/proteomics_module/fragpipe-run.py",
                                 " --db_path ", shQuote(renamed_fragpipe_prot_db_file),
                                 " --mass_spec_info_path ", shQuote(mass_spec_info_file),
                                 " --output_dir ", shQuote(output_dir),
                                 " --protease1 ", shQuote(protease1))

  if (protease2 != "none") {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --protease2 ", shQuote(protease2))
  }

  if (input$user_add_contaminants) {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --add_contaminants")
  }

  if (input$user_perform_quantification) {
    command_run_fragpipe <- paste0(command_run_fragpipe, " --perform_quantification")
  }

  command_run_fragpipe <- paste0(command_run_fragpipe,
                                 " --fragpipe_path ", shQuote("/home/user/Desktop/GenomeProt/fragpipe-23.1/"),
                                 " --num_threads ", shQuote(fragpipe_cpu_threads),
                                 " --memory_limit ", shQuote(fragpipe_memory_limit))
  message(command_run_fragpipe)
  system(command_run_fragpipe)

  # Find files to zip
  files_to_zip_fragpipe <- c("peptide.tsv")
  if (input$user_perform_quantification) {
    files_to_zip_fragpipe <- c(files_to_zip_fragpipe, file.path("dia-quant-output", "report.pr_matrix.tsv"))
  }

  error_msg <- ""
  if (base::all(file.exists(files_to_zip_fragpipe))) {
    zipfile_path_fragpipe <- file.path("..", "fragpipe_results.zip")
    zip(zipfile = zipfile_path_fragpipe, files = files_to_zip_fragpipe)
  } else {
    error_msg <- "Error: Proteomics search failed. Please check the proteomics module log file generated in the output directory for details."
  }

  # Go back to the old directory
  setwd(old_wd)

  return(error_msg)
}

is_peptide_reformat_input_valid <- function(input) {
  error_msg <- ""

  # Check if the user has uploaded proteomics results files yet
  user_orig_proteomics_files <- input$user_orig_proteomics_files
  if (is.null(user_orig_proteomics_files)) {
    error_msg <- "Error: At least one proteomics results file must be uploaded."
    return(error_msg)
  }

  # Every proteomics results file must have a file extension of '.txt', '.csv' or '.tsv' and must not be empty
  num_files <- nrow(user_orig_proteomics_files)
  filenames <- user_orig_proteomics_files$name
  for (i in 1:num_files) {
    filename <- filenames[i]
    if (!any(endsWith(filename, c(".txt", ".csv", ".tsv")))) {
      error_msg <- "Error: Please ensure that all uploaded proteomics results files have a file extension of '.txt', '.csv' or '.tsv'."
      return(error_msg)
    }

    file_path <- user_orig_proteomics_files$datapath[i]
    if (file.size(file_path) == 0) {
      error_msg <- "Error: Please ensure that none of the uploaded proteomics results files are empty."
      return(error_msg)
    }
  }

  # The proteomics search tool selected must be one of 'Spectronaut', 'FragPipe' or 'FragPipe_quant'
  search_tool <- input$proteomics_search_tool
  if (is.null(search_tool) || !(search_tool %in% c("Spectronaut", "FragPipe", "FragPipe_quant"))) {
    error_msg <- "Error: Please select one of the proteomics search tools listed."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

integration_reformat_server <- function(input, session) {
  session_id <- session$token                                               # session ID
  outdir_peptide_results <- file.path(session_id, "peptide_results_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_peptide_results)) {
    dir.create(outdir_peptide_results)
  }

  user_orig_proteomics_files <- input$user_orig_proteomics_files
  search_tool <- input$proteomics_search_tool
  orig_folder <- dirname(user_orig_proteomics_files$datapath[1])

  # move all uploaded proteomics files into the output directory
  num_files <- nrow(user_orig_proteomics_files)
  new_paths <- c()
  for (i in 1:num_files) {
    orig_name <- user_orig_proteomics_files$name[i]

    new_name <- stringr::str_replace(orig_name, "\\.txt$", ".tsv")
    new_name <- stringr::str_replace(orig_name, "\\.csv$", ".tsv")
    if (new_name == "peptide_data.tsv") {
      new_name <- "peptide_data_renamed.tsv"
    }

    new_path <- file.path(orig_folder, new_name)
    new_paths <- c(new_paths, new_path)
  }

  mapply(file.rename, from = user_orig_proteomics_files$datapath, to = new_paths)

  command_peptide_reformat <- paste0(conda_command, "Rscript bin/integration_module/reformat_peptide_results.R",
                                     " -d ", shQuote(orig_folder),
                                     " -s ", shQuote(search_tool))

  # reformat the proteomics results
  message(command_peptide_reformat)
  system(command_peptide_reformat)

  error_msg <- ""
  reformatted_file <- file.path(orig_folder, "peptide_data.tsv")
  if (!file.exists(reformatted_file)) {
    error_msg <- "Error: Failed to reformat proteomics results. Please ensure that the files are formatted correctly and the correct search tool is provided."
  } else {
    output_location <- file.path(outdir_peptide_results, "peptide_data.tsv")
    file.copy(reformatted_file, output_location)
  }

  return(error_msg)
}

is_integration_input_valid <- function(input) {
  error_msg <- ""

  # Check if the user has uploaded a reformatted proteomics results file yet
  user_proteomics_file <- input$user_proteomics_file
  if (is.null(user_proteomics_file)) {
    error_msg <- "Error: Please upload a reformatted proteomics results file to perform proteogenomics integration."
    return(error_msg)
  }

  # The reformatted proteomics results file must have a file extension of '.txt', '.csv' or '.tsv'
  filename <- user_proteomics_file$name
  if (!any(endsWith(filename, c(".txt", ".csv", ".tsv")))) {
    error_msg <- "Error: The reformatted proteomics results file must have a file extension of '.txt', '.csv' or '.tsv'."
    return(error_msg)
  }

  # The reformatted proteomics results file must not be empty
  file_path <- user_proteomics_file$datapath
  if (file.size(file_path) == 0) {
    error_msg <- "Error: The reformatted proteomics results file must not be empty."
    return(error_msg)
  }

  # Check if the user has uploaded proteome_database_metadata.txt yet
  user_metadata_file <- input$user_metadata_file
  if (is.null(user_metadata_file)) {
    error_msg <- "Error: Please upload 'proteome_database_metadata.txt' to perform proteogenomics integration."
    return(error_msg)
  }

  # proteome_database_metadata.txt must have a file extension of '.txt'
  filename <- user_metadata_file$name
  if (!endsWith(filename, ".txt")) {
    error_msg <- "Error: 'proteome_database_metadata.txt' must have a file extension of '.txt'."
    return(error_msg)
  }

  # proteome_database_metadata.txt must not be empty
  file_path <- user_metadata_file$datapath
  if (file.size(file_path) == 0) {
    error_msg <- "Error: proteome_database_metadata.txt must not be empty."
    return(error_msg)
  }

  # Check if the user has uploaded proteome_database_transcripts.gtf yet
  user_post_gtf_file <- input$user_post_gtf_file
  if (is.null(user_post_gtf_file)) {
    error_msg <- "Error: Please upload 'proteome_database_transcripts.gtf' to perform proteogenomics integration."
    return(error_msg)
  }

  # proteome_database_transcripts.gtf must have a file extension of '.gtf'
  filename <- user_post_gtf_file$name
  if (!endsWith(filename, ".gtf")) {
    error_msg <- "Error: 'proteome_database_transcripts.gtf' must have a file extension of '.gtf'."
    return(error_msg)
  }

  # proteome_database_transcripts.gtf must not be empty
  file_path <- user_post_gtf_file$datapath
  if (file.size(file_path) == 0) {
    error_msg <- "Error: proteome_database_transcripts.gtf must not be empty."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

integration_server <- function(input, session) {
  session_id <- session$token                           # session ID
  outdir_integ <- file.path(session_id, "integ_output") # output directory

  # create the output directory
  if (!dir.exists(outdir_integ)) {
    dir.create(outdir_integ)
  }

  command_map_peptides <- paste0(conda_command, "Rscript bin/integration_module/map_peptides_generate_outputs.R",
                                 " -p ", shQuote(input$user_proteomics_file$datapath),
                                 " -m ", shQuote(input$user_metadata_file$datapath),
                                 " -g ", shQuote(input$user_post_gtf_file$datapath),
                                 " -s ", shQuote(outdir_integ))

  # map peptides and perform integration
  message(command_map_peptides)
  system(command_map_peptides)

  error_msg <- ""
  expected_files <- file.path(outdir_integ, c("peptides.bed12", "ORFs.bed12", "transcripts.bed12", "transcripts_and_ORFs_for_isovis.gtf", "peptide_info.tsv", "combined_annotations.gtf"))
  if (!base::all(file.exists(expected_files))) {
    error_msg <- "Error: Proteogenomics integration failed. Some results files are missing. Please ensure the input files you have uploaded are valid."
    return(error_msg)
  }

  top_level_dir <- getwd()
  report_rmd_file <- file.path(top_level_dir, "bin", "integration_module", "integration_summary_report.Rmd")
  report_outdir <- file.path(top_level_dir, outdir_integ)

  command_create_report <- paste0(conda_command, "Rscript bin/integration_module/generate_integration_summary_report.R",
                                  " -i ", shQuote(report_rmd_file),
                                  " -o ", shQuote(report_outdir))

  # create the summary report
  message(command_create_report)
  system(command_create_report)

  # create a zip file with results
  files_to_zip_int <- c("summary_report.html", "peptide_info.tsv", "report_images/",
                        "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                        "peptides.bed12", "ORFs.bed12", "transcripts.bed12", "ncORF_stats.xlsx")

  # set the path to the ZIP file (in the session_id directory)
  zipfile_path_int <- file.path("..", "integration_results.zip")

  # temp change the working dir to outdir_integ
  setwd(outdir_integ)

  if (base::all(file.exists(files_to_zip_int))) {
    if (file.exists(zipfile_path_int)) {
      file.remove(zipfile_path_int)
    }

    # zip files
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
  } else {
    error_msg <- "Error: Failed to generate proteogenomics integration summary report. Please ensure the input files you have uploaded are valid."
  }

  # go back to starting dir
  setwd(top_level_dir)

  return(error_msg)
}

# main shiny app server
server <- function(input, output, session) {
  session_id <- session$token   # session ID
  message(paste0("Session: ", session_id))

  # create session id tmp directory each time app is run
  if (!dir.exists(session_id)) {
    dir.create(session_id)
  }

  # TEST DATA

  # create reactive value for the test data results file
  file_available_test_data_output <- reactiveVal(FALSE)

  # run test data through GenomeProt when the submit button is pressed
  observeEvent(input$test_data_submit_button, {
    file_available_test_data_output(FALSE)

    message("Validating input for running test data...")
    error_msg <- is_test_data_input_valid()
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "test-data-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "test-data-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "test_data_submit_button", spinnerId = "test-data-loading-container")) # disable submit button

    # run the test data through GenomeProt
    test_data_output_file <- file.path(session_id, "test_data_output", "test_data_results.zip")
    if (file.exists(test_data_output_file)) {
      file.remove(test_data_output_file)
    }
    error_msg <- test_data_server(session)

    session$sendCustomMessage("enableButton", list(id = "test_data_submit_button", spinnerId = "test-data-loading-container"))

    # check if the test data results file is created
    if (!file.exists(test_data_output_file) & (error_msg == "")) {
      error_msg <- "Error: Failed to run test data through GenomeProt. Please ensure that you have followed the installation instructions and that there are sufficient system resources to run the test data through."
    }

    if (error_msg != "") {
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "test-data-status-msg-container", color = "red"))
    } else {
      file_available_test_data_output(TRUE)
    }
  })

  # enable download once files are available
  observe({
    if (file_available_test_data_output()) {
      shinyjs::enable("test_data_download_button")
      shinyjs::runjs("document.getElementById('test_data_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "test_data_submit_button", spinnerId = "test-data-loading-container")) # re-enable submit button
    }
  })

  # download handler for the test data results file
  output$test_data_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_test_data_results.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "test_data_output", "test_data_results.zip"), file)
    }
  )

  # DATABASE MODULE

  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)

  # generate the database when the submit button is pressed
  observeEvent(input$db_submit_button, {
    file_available_db(FALSE)

    message("Validating database generation module input...")
    error_msg <- is_database_generation_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "db-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "db-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "db_submit_button", spinnerId = "db-loading-container")) # disable submit button

    # run the database generation server
    results_zip_file <- file.path(session_id, "database_results.zip")
    if (file.exists(results_zip_file)) {
      file.remove(results_zip_file)
    }

    # TODO: Add error messages for the database generation module servers
    if (input$input_type == "fastq_input" & input$sequencing_type == "long-read") {
      fastq_server(input, session)
      bambu_server(input, session)
      database_server(input, session)
    } else if (input$input_type == "fastq_input" & input$sequencing_type == "short-read") {
      fastq_server(input, session)
      database_server(input, session)
    } else if (input$input_type == "bam_input" & input$sequencing_type == "long-read") {
      bambu_server(input, session)
      database_server(input, session)
    } else if (input$input_type == "bam_input" & input$sequencing_type == "short-read") {
      bam_server(input, session)
      database_server(input, session)
    } else {
      database_server(input, session)
    }

    session$sendCustomMessage("enableButton", list(id = "db_submit_button", spinnerId = "db-loading-container"))

    # check if the zip file is created
    error_msg <- ""
    if (!file.exists(results_zip_file) & (error_msg == "")) {
      error_msg <- "Error: Database generation failed. Please ensure the files you have uploaded are formatted correctly and you have selected the correct input type and sequencing type."
    }

    if (error_msg != "") {
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "db-status-msg-container", color = "red"))
    } else {
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
                               paste0("Upload mass spectrometry file #", mass_spec_file_num, ":"),
                               accept = c(".mzML", ".mzXML", ".mgf", ".mzBIN", ".raw", ".d"))
              ),
              column(4,
                     selectInput(paste0("mass_spec_file_data_type_input_", mass_spec_file_num),
                                 "Data type:",
                                 choices = c("DDA" = "DDA",
                                             "DDA+" = "DDA+",
                                             "DIA" = "DIA",
                                             "DIA-Quant" = "DIA-Quant",
                                             "DIA-Lib" = "DIA-Lib",
                                             "GPF-DIA" = "GPF-DIA"))
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
    results_zip_file <- file.path(session_id, "fragpipe_results.zip")
    if (file.exists(results_zip_file)) {
      file.remove(results_zip_file)
    }
    error_msg <- fragpipe_server(input, session_id, mass_spec_file_num)

    session$sendCustomMessage("enableButton", list(id = "fragpipe_submit_button", spinnerId = "fragpipe-loading-container"))

    # check if the zip file is created
    if (!file.exists(results_zip_file) & (error_msg == "")) {
      error_msg <- "Error: Proteomics search failed. Please check the proteomics module log file generated in the output directory for details."
    }

    if (error_msg != "") {
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

  # create reactive value for the reformatted proteomics results file
  file_available_integ_reformatted <- reactiveVal(FALSE)

  # run integration function when submit is pressed
  observeEvent(input$integ_reformat_submit_button, {
    file_available_integ_reformatted(FALSE)

    message("Validating input for reformatting proteomics results...")
    error_msg <- is_peptide_reformat_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-reformat-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "integ-reformat-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "integ_reformat_submit_button", spinnerId = "integ-reformat-loading-container")) # disable submit button

    # run the proteomics results reformatting server
    reformatted_file <- file.path(session_id, "peptide_results_output", "peptide_data.tsv")
    if (file.exists(reformatted_file)) {
      file.remove(reformatted_file)
    }
    error_msg <- integration_reformat_server(input, session)

    session$sendCustomMessage("enableButton", list(id = "integ_reformat_submit_button", spinnerId = "integ-reformat-loading-container"))

    # check if the reformatted proteomics results TSV file is created
    if (!file.exists(reformatted_file) & (error_msg == "")) {
      error_msg <- "Error: Failed to reformat proteomics results. Please ensure that the files are formatted correctly and the correct search tool is provided."
    }

    if (error_msg != "") {
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-reformat-status-msg-container", color = "red"))
    } else {
      file_available_integ_reformatted(TRUE)
    }
  })

  # enable download once files are available
  observe({
    if (file_available_integ_reformatted()) {
      shinyjs::enable("integ_reformat_download_button")
      shinyjs::runjs("document.getElementById('integ_reformat_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_reformat_submit_button", spinnerId = "integ-reformat-loading-container")) # re-enable submit button
    }
  })

  # download handler for the reformatted proteomics results TSV file
  output$integ_reformat_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_peptide_data.tsv")
    },
    content = function(file) {
      file.copy(file.path(session_id, "peptide_results_output", "peptide_data.tsv"), file)
    }
  )

  # create reactive value for the database zip
  file_available_integ <- reactiveVal(FALSE)

  # run integration function when submit is pressed
  observeEvent(input$integ_submit_button, {
    file_available_integ(FALSE)

    message("Validating integration module input...")
    error_msg <- is_integration_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "integ-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button

    # run integration server
    integration_results_zip <- file.path(session_id, "integration_results.zip")
    if (file.exists(integration_results_zip)) {
      file.remove(integration_results_zip)
    }
    error_msg <- integration_server(input, session)

    session$sendCustomMessage("enableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container"))

    # check if the zip file is created
    if (!file.exists(integration_results_zip) & (error_msg == "")) {
      error_msg <- "Proteogenomics integration failed. Some results files are missing. Please ensure the input files you have uploaded are valid."
    }

    if (error_msg != "") {
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-status-msg-container", color = "red"))
    } else {
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

  # download handler for the integration module results zip
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