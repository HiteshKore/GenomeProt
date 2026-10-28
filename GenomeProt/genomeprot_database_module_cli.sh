#!/bin/bash

# Get the path to the conda binary
conda_path=$(/usr/bin/which conda)
if [ -z "$conda_path" ]; then
    echo 'Error: Failed to find the conda binary. Check if conda is installed or if it is within a directory specified in the $PATH variable.' >&2
    exit 1
fi

# Name of the Conda environment
genomeprot_env="GenomeProt_env"

# Check if the Conda environment exists
if ! ("$conda_path" env list | grep -q "^$genomeprot_env\s"); then
    echo "Error: The Conda environment '$genomeprot_env' does not exist." >&2
    exit 1
fi

function process_short_read_samples() {
    fq_r1_pe=("${!1}")
    fq_r2_pe=("${!2}")
    fq_se=("${!3}")
    bam_files=("${!4}")
    genome_fa=$5
    transcriptome_fa=$6
    ref_gtf=$7
    output_dir=$8
    threads=$9

    # If paired-end and/or single-end FASTQ samples are provided, generate a salmon index before processing them
    if [[ ${#fq_r1_pe[@]} -gt 0 ]] || [[ ${#fq_se[@]} -gt 0 ]]; then
        echo 'Generating salmon index' | tee -a "$log_file"

        decoy_file="${output_dir}decoys.txt"
        salmon_index_dir="${output_dir}salmon_index"

        # If both the reference genome and reference transcriptome are gzipped, the 'gentrome' (concatenation of
        # the genome and the transcriptome) can be gzipped
        if [[ "$transcriptome_fa" == *.gz ]] && [[ "$genome_fa" == *.gz ]]; then
            echo 'Processing gzipped files' | tee -a "$log_file"

            # Generate decoy file (sequence names without '>')
            grep '^>' <(gzip -c -d "$genome_fa") | cut -d ' ' -f 1 | sed 's/>//g' > "$decoy_file"

            # Concatenate reference and transcriptome to create gentrome
            cat "$transcriptome_fa" "$genome_fa" > "${output_dir}gentrome.fa.gz"

            # Run Salmon index
            "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                -t "${output_dir}gentrome.fa.gz" -d "$decoy_file" \
                -p "$threads" -i "$salmon_index_dir"

        # If both the reference genome and reference transcriptome are uncompressed, the gentrome is uncompressed
        elif [[ "$transcriptome_fa" != *.gz ]] && [[ "$genome_fa" != *.gz ]]; then
            echo 'Processing uncompressed files' | tee -a "$log_file"

            # Generate decoy file (sequence names without '>')
            grep '^>' "$genome_fa" | cut -d ' ' -f 1 | sed 's/>//g' > "$decoy_file"

            # Concatenate reference and transcriptome to create gentrome
            cat "$transcriptome_fa" "$genome_fa" > "${output_dir}gentrome.fa"

            # Run Salmon index
            "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                -t "${output_dir}gentrome.fa" -d "$decoy_file" \
                -p "$threads" -i "$salmon_index_dir"

        # If only the reference transcriptome is gzipped, decompress it first
        elif [[ "$transcriptome_fa" == *.gz ]] && [[ "$genome_fa" != *.gz ]]; then
            echo 'Processing files' | tee -a "$log_file"

            # Generate decoy file (sequence names without '>')
            grep '^>' "$genome_fa" | cut -d ' ' -f 1 | sed 's/>//g' > "$decoy_file"

            # Create a temporary decompressed reference transcriptome
            temp_transcriptome_fa="${output_dir}decompressed_transcriptome.fa"
            gzip -c -d "$transcriptome_fa" > "$temp_transcriptome_fa"

            # Concatenate reference and transcriptome to create gentrome
            cat "$temp_transcriptome_fa" "$genome_fa" > "${output_dir}gentrome.fa"

            # Run Salmon index
            "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                -t "${output_dir}gentrome.fa" -d "$decoy_file" \
                -p "$threads" -i "$salmon_index_dir"

            # Delete the temporary decompressed reference transcriptome
            rm "$temp_transcriptome_fa"

        # If only the reference genome is gzipped, decompress it first
        else
            echo 'Processing files' | tee -a "$log_file"

            # Create a temporary decompressed reference genome
            temp_genome_fa="${output_dir}decompressed_genome.fa"
            gzip -c -d "$genome_fa" > "$temp_genome_fa"

            # Generate decoy file (sequence names without '>')
            grep '^>' "$temp_genome_fa" | cut -d ' ' -f 1 | sed 's/>//g' > "$decoy_file"

            # Concatenate reference and transcriptome to create gentrome
            cat "$transcriptome_fa" "$temp_genome_fa" > "${output_dir}gentrome.fa"

            # Run Salmon index
            "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                -t "${output_dir}gentrome.fa" -d "$decoy_file" \
                -p "$threads" -i "$salmon_index_dir"

            # Note: Do NOT delete the temporary decompressed reference genome as it will be used for generating the proteome database later
        fi

        # Process paired-end samples
        if [[ ${#fq_r1_pe[@]} -gt 0 ]]; then
            for paired_file_index in $(seq 0 $(( ${#fq_r1_pe[@]} - 1 ))); do
                r1_file_path=${fq_r1_pe[$paired_file_index]}
                r2_file_path=${fq_r2_pe[$paired_file_index]}
                sample_id=$(basename "$r1_file_path" | sed -E 's/_R1\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')

                echo "Quantifying paired-end reads from $sample_id" | tee -a "$log_file"

                "$conda_path" run -n "$genomeprot_env" salmon quant \
                    -i "$salmon_index_dir" -p "$threads" -l A \
                    -1 "$r1_file_path" -2 "$r2_file_path" \
                    -o "${output_dir}${sample_id}"
            done

            echo 'Paired-end samples processed.' | tee -a "$log_file"
        fi

        # Process single-end samples
        if [[ ${#fq_se[@]} -gt 0 ]]; then
            for file_path in "${fq_se[@]}"; do
                sample_id=$(basename "$file_path" | sed -E 's/\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
                echo "Quantifying single-end reads from $sample_id" | tee -a "$log_file"

                "$conda_path" run -n "$genomeprot_env" salmon quant \
                    -i "$salmon_index_dir" -p "$threads" -l A \
                    -r "$file_path" \
                    -o "${output_dir}${sample_id}"
            done

            echo 'Single-end samples processed.' | tee -a "$log_file"
        fi

    # Otherwise, process BAM samples
    else
        genome_file="$genome_fa"

        # If the reference genome is gzipped, decompress it first
        if [[ "$genome_fa" == *.gz ]]; then
            genome_file="${output_dir}decompressed_genome.fa"
            gzip -c -d "$genome_fa" > "$genome_file"
        fi

        # Generate the transcriptome file
        transcriptome_file="${output_dir}transcript.fa"
        "$conda_path" run -n "$genomeprot_env" gffread -w "$transcriptome_file" -g "$genome_file" "$ref_gtf"

        # Note: Do NOT delete the temporary decompressed reference genome as it will be used for generating the proteome database later

        for file_path in "${bam_files[@]}"; do
            sample_id=$(basename "$file_path" | sed -E 's/\.bam$//g')
            echo "Quantifying reads from $sample_id" | tee -a "$log_file"

            "$conda_path" run -n "$genomeprot_env" salmon quant \
                -t "$transcriptome_file" -p "$threads" -l A \
                -a "$file_path" \
                -o "${output_dir}${sample_id}"
        done

        echo 'BAM samples processed.' | tee -a "$log_file"
    fi

    echo 'Generating transcript counts from Salmon quantification files' | tee -a "$log_file"

    "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/matrix_compilation_salmon.R -s "$output_dir" -g "$ref_gtf" -o "${output_dir}salmon_transcript_counts.txt"

    echo 'Transcript counts generated.' | tee -a "$log_file"
}

function process_long_read_samples() {
    fastq_files=("${!1}")
    genome_fa=$2
    ref_gtf=$3
    sample_dir=$4
    output_dir=$5
    threads=$6
    organism=$7

    bam_folder="$sample_dir"

    # If FASTQ samples are provided, map the reads to the reference genome and generate BAM files
    if [[ ${#fastq_files[@]} -gt 0 ]]; then
        genome_index_filename_without_extension=$(basename "$genome_fa" | sed -E 's/\.(fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
        genome_index_path="${output_dir}/${genome_index_filename_without_extension}.mmi"

        echo 'Generating minimap2 index' | tee -a "$log_file"
        "$conda_path" run -n "$genomeprot_env" minimap2 -ax splice:hq -d "$genome_index_path" "$genome_fa"
        echo 'Indexing complete.' | tee -a "$log_file"

        bam_folder="$output_dir"

        echo 'Genome alignment started' | tee -a "$log_file"
        for fastq_file_path in "${fastq_files[@]}"; do
            file_prefix=$(basename "$fastq_file_path" | sed -E 's/\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
            echo "Processing: $file_prefix" | tee -a "$log_file"
            "$conda_path" run -n "$genomeprot_env" bash -c "minimap2 --sam-hit-only --secondary=no -t $threads -ax splice:hq $genome_index_path $fastq_file_path | samtools view -bh | samtools sort -@ $threads -o ${bam_folder}${file_prefix}.bam"
        done
        echo 'Genome alignment completed.' | tee -a "$log_file"
    fi

    echo 'Performing transcript identification and quantification with Bambu' | tee -a "$log_file"
    "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/run_bambu.R -b "$bam_folder" -g "$ref_gtf" -t "$threads" -s "$organism" -o "$output_dir" 2>&1 | tee -a "$log_file"
    echo 'Transcripts identified and transcript counts generated.' | tee -a "$log_file"

    custom_gtf_file="${output_dir}bambu_transcript_annotations.gtf"
    tx_counts_file="${output_dir}bambu_transcript_counts.txt"

    # Rename files
    mv "${output_dir}extended_annotations.gtf" "$custom_gtf_file"
    mv "${output_dir}counts_transcript.txt" "$tx_counts_file"
    mv "${output_dir}counts_gene.txt" "${output_dir}bambu_gene_counts.txt"

    echo 'Running gffcompare on the reference transcriptome GTF and the Bambu transcript GTF annotation' | tee -a "$log_file"
    "$conda_path" run -n "$genomeprot_env" gffcompare -r "$ref_gtf" "$custom_gtf_file"
    mv "${output_dir}/gffcmp.bambu_transcript_annotations.gtf.tmap" "${output_dir}/gffcompare.tmap.txt"
    rm gffcmp*
    rm "${output_dir}"/gffcmp*
    echo 'gffcompare finished running.' | tee -a "$log_file"
}

function generate_database() {
    custom_gtf=$1
    output_dir=$2
    ref_gtf=$3
    organism=$4
    min_orf_length=$5
    vcf_fn=$6
    orf_input_type=$7
    five_utr_orf=$8
    three_utr_orf=$9
    ref_proteome=${10}
    ref_genome=${11}
    tx_counts_file=${12}
    min_tx_count=${13}
    input_type=${14}
    threads=${15}

    # No VCF file was provided
    if [ "$vcf_fn" = "None" ]; then
        echo 'Generating the proteome database' | tee -a "$log_file"

        if [ "$tx_counts_file" = "None" ]; then
            echo 'No transcript counts file provided.' | tee -a "$log_file"
            "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/generate_proteome.R -g "$custom_gtf" -r "$ref_gtf" -o "$organism" -l "$min_orf_length" -u "$five_utr_orf" -d "$three_utr_orf" -s "$output_dir"
        else
            echo "Transcript counts file: $tx_counts_file" | tee -a "$log_file"
            echo "Transcript expression cutoff: $min_tx_count" | tee -a "$log_file"
            "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/generate_proteome.R -g "$custom_gtf" -r "$ref_gtf" -c "$tx_counts_file" -m "$min_tx_count" -o "$organism" -l "$min_orf_length" -u "$five_utr_orf" -d "$three_utr_orf" -s "$output_dir"
        fi

        echo 'Annotating the proteome database' | tee -a "$log_file"

        if [ "$input_type" = "GTF" ]; then
            "$conda_path" run -n "$genomeprot_env" python ./bin/database_module/annotate_proteome.py "$ref_gtf" "$ref_proteome" "${output_dir}ORFome_aa.txt" "$custom_gtf" "$output_dir" "$orf_input_type" "$min_orf_length" None "$organism" "$threads" 2000
        else
            "$conda_path" run -n "$genomeprot_env" python ./bin/database_module/annotate_proteome.py "$ref_gtf" "$ref_proteome" "${output_dir}ORFome_aa.txt" "${output_dir}proteome_database_transcripts.gtf" "$output_dir" "$orf_input_type" "$min_orf_length" None "$organism" "$threads" 2000
        fi

    # A VCF file was provided
    else
        echo 'Generating the variant-aware proteome database' | tee -a "$log_file"

        if [ "$tx_counts_file" = "None" ]; then
            echo 'No transcript counts file provided.' | tee -a "$log_file"
            "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/generate_proteome.R -G "$ref_genome" -g "$custom_gtf" -r "$ref_gtf" -o "$organism" -l "$min_orf_length" -u "$five_utr_orf" -d "$three_utr_orf" -v "$vcf_fn" -s "$output_dir"
        else
            echo "Transcript counts file: $tx_counts_file" | tee -a "$log_file"
            echo "Transcript expression cutoff: $min_tx_count" | tee -a "$log_file"
            "$conda_path" run -n "$genomeprot_env" Rscript ./bin/database_module/generate_proteome.R -G "$ref_genome" -g "$custom_gtf" -r "$ref_gtf" -c "$tx_counts_file" -m "$min_tx_count" -o "$organism" -l "$min_orf_length" -u "$five_utr_orf" -d "$three_utr_orf" -v "$vcf_fn" -s "$output_dir"
        fi

        echo 'Annotating the variant-aware proteome database' | tee -a "$log_file"

        mutant_db="${output_dir}Mutant_ORFome_aa.txt"
        "$conda_path" run -n "$genomeprot_env" python ./bin/database_module/annotate_proteome.py "$ref_gtf" "$ref_proteome" "${output_dir}ORFome_aa.txt" "${output_dir}proteome_database_transcripts.gtf" "$output_dir" "$orf_input_type" "$min_orf_length" "$mutant_db" "$organism" "$threads" 2000
    fi

    echo 'Proteome database generated and annotated.' | tee -a "$log_file"
}

function is_uint() {
    case "$1" in
        ''|*[!0-9]*) return 1 ;;
    esac
}

function usage() {
    cat << EOF
Program: GenomeProt proteome database generation module (command-line)
Version: 0.0.1v

 ▗▄▄▖▗▄▄▄▖▗▖  ▗▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▖ ▗▄▄▖  ▗▄▖▗▄▄▄▖
▐▌   ▐▌   ▐▛▚▖▐▌▐▌ ▐▌▐▛▚▞▜▌▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌ █
▐▌▝▜▌▐▛▀▀▘▐▌ ▝▜▌▐▌ ▐▌▐▌  ▐▌▐▛▀▀▘▐▛▀▘ ▐▛▀▚▖▐▌ ▐▌ █
▝▚▄▞▘▐▙▄▄▖▐▌  ▐▌▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▌   ▐▌ ▐▌▝▚▄▞▘ █

Usage: genomeprot_database_module_cli.sh [-h] [-p <threads>] [-O <output directory>] [-s <sequencing type>] [-t <data type>] [-d <sample directory>] [-o <organism>] [-g <reference genome FASTA>] [-r <reference transcriptome FASTA>] [-a <reference transcriptome GTF>] [-G <custom transcript annotation GTF>] [-v <VCF file>] [-T <ORF type>] [-l <minimum ORF length (AA)>] [-U <find upstream / 5' UTR ORFs?>] [-D <find downstream / 3' UTR ORFs?>] [-m <minimum transcript count>] [-C <transcript counts file>]

Generate a proteome database from long-read or short-read RNA-seq data.
Supported input data types: FASTQ, BAM or GTF.

GENERAL
  -h, --help                      Display this help and exit
  -p, --threads <int>             Number of threads (default: 1)
  -O, --output-dir <dir>          Output directory [required]

INPUT
  -s, --sequencing-type <type>    Sequencing type: long-read | short-read [required]
  -t, --data-type <type>          Input data type: FASTQ | BAM | GTF [required]
  -d, --sample-dir <dir>          Sample directory [required for FASTQ and BAM input; FASTQ files can be gzipped]
  -o, --organism <name>           Organism: HUMAN | CAEEL | DROME | MOUSE | RAT | DANRE [required]

REFERENCES
  -g, --genome <file>             Reference genome FASTA [required for FASTQ and BAM input and for generating a variant-aware proteome database; can be gzipped]
  -r, --transcriptome <file>      Reference transcriptome FASTA [required for short-read RNA-seq FASTQ input; can be gzipped]
  -a, --reference-gtf <file>      Reference transcriptome GTF [required; must contain GENCODE annotations only]
  -G, --custom-gtf <file>         Custom transcript annotation GTF [required for long-read GTF input]
  -v, --vcf <file>                VCF file with variants from all samples (optional; used for generating a variant-aware proteome database)

ORF OPTIONS
  -T, --orf-type <type>           ORF type: canonical | all (default: all)
  -l, --min-orf-length <int>      Minimum ORF length in amino acids (default: 30)
  -U, --upstream-orf <bool>       Include upstream (5' UTR) ORFs: TRUE | FALSE (default: FALSE)
  -D, --downstream-orf <bool>     Include downstream (3' UTR) ORFs: TRUE | FALSE (default: FALSE)

EXPRESSION FILTERING
  -m, --min-tx-count <int>        Minimum transcript expression cutoff (default: 5)
  -C, --tx-counts <file>          Transcript counts file [required for short-read GTF input]

=== Expected FASTQ and BAM sample formats ===

This script searches for FASTQ and BAM samples in the sample directory specified. As the script does not perform recursive search, samples placed in folders within that directory are ignored.
For example, if the FASTQ sample directory specified contains 'sample1.fastq', 'sample2.fq.gz', 'folder1/sample3.fastq' and 'folder1/folder2/sample4.fastq', only the first two samples will be used.

FASTQ samples must have a file extension of '.fastq', '.fq', '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped.
Long-read and single-end short-read FASTQ samples must follow the format '<sample ID><file extension>'.
Paired-end short-read FASTQ samples must follow the format '<sample ID>_R1<file extension>' or '<sample ID>_R2<file extension>', with the same sample ID provided for both ends (R1 and R2).

Sample IDs must not be empty. Any files with empty sample IDs are ignored (e.g. '.fastq.gz', '_R1.fastq', '_R2.fa').
Sample IDs must be unique. No duplicate sample IDs are allowed for the following cases:
1. All long-read FASTQ samples
2. Either end of the paired-end short-read FASTQ samples (R1 or R2)
3. Single-end short-read FASTQ samples
4. All short-read FASTQ samples (including paired-end and single-end)

BAM samples must have a file extension of '.bam'.

=== Expected file formats ===

The reference genome FASTA (-g / --genome) and reference transcriptome FASTA (-r / --transcriptome) must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped.
The reference transcriptome GTF (-a / --reference-gtf) and custom transcript annotation GTF (-G / --custom-gtf) must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'.
The VCF file (-v / --vcf) must have a file extension of '.vcf'.
The transcript counts file (-C / --tx-counts) must have a file extension of '.txt', '.csv' or '.tsv'.

EOF
}

# Ensure that the script exits if an error occurs
set -e

# Check if no arguments were provided
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Parse options using getopt
OPTS=$(getopt -o s:o:l:t:g:p:a:G:v:T:r:U:D:m:C:d:O:h \
              -l sequencing-type:,organism:,min-orf-length:,data-type:,genome:,threads:,reference-gtf:,custom-gtf:,vcf:,orf-type:,transcriptome:,upstream-orf:,downstream-orf:,min-tx-count:,tx-counts:,sample-dir:,output-dir:,help \
              -n "genomeprot_database_module_cli.sh" -- "$@")

if [ $? -ne 0 ]; then
    echo 'Error: Failed to parse options.' >&2
    usage
    exit 1
fi

# Re-evaluate the set of positional parameters
eval set -- "$OPTS"

# Process options
while true; do
    case "$1" in
        -h | --help )               usage; exit 0 ;;
        -p | --threads )            threads="$2"; shift 2 ;;
        -O | --output-dir )         output_directory="$2"; shift 2 ;;
        -s | --sequencing-type )    sequencing_platform="$2"; shift 2 ;;
        -t | --data-type )          data_type="$2"; shift 2 ;;
        -d | --sample-dir )         sample_directory="$2"; shift 2 ;;
        -o | --organism )           organism="$2"; shift 2 ;;
        -g | --genome )             genome_fa="$2"; shift 2 ;;
        -r | --transcriptome )      transcriptome_fa="$2"; shift 2 ;;
        -a | --reference-gtf )      reference_gtf="$2"; shift 2 ;;
        -G | --custom-gtf )         custom_gtf="$2"; shift 2 ;;
        -v | --vcf )                vcf="$2"; shift 2 ;;
        -T | --orf-type )           orf_type="$2"; shift 2 ;;
        -l | --min-orf-length )     min_orf_length="$2"; shift 2 ;;
        -U | --upstream-orf )       upstream_orf="$2"; shift 2 ;;
        -D | --downstream-orf )     downstream_orf="$2"; shift 2 ;;
        -m | --min-tx-count )       min_tx_count="$2"; shift 2 ;;
        -C | --tx-counts )          tx_counts_file="$2"; shift 2 ;;
        -- ) shift; break ;;        # End of options
        * ) echo "Invalid option: $1"; usage; exit 1 ;;
    esac
done

# Check whether the input arguments are correct

if [ -z "$output_directory" ]; then
    echo 'Error: An output directory must be provided.' >&2
    exit 1
fi

if [[ ! "$output_directory" =~ /$ ]]; then
    output_directory="$output_directory/"
fi

# If the output directory does not exist yet, create it
if ! [ -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

log_file="${output_directory}database_module_logfile.txt"
echo "Log file: $log_file" | tee -a "$log_file"

if [ -z "$sequencing_platform" ]; then
    echo 'Error: A sequencing type must be specified.' | tee -a "$log_file" >&2
    exit 1
elif [ "$sequencing_platform" = "long-read" ] || [ "$sequencing_platform" = "short-read" ]; then
    echo "Sequencing type: $sequencing_platform" | tee -a "$log_file"
else
    echo "Error: The sequencing type '$sequencing_platform' is not supported. Please specify either 'long-read' or 'short-read'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$data_type" ]; then
    echo 'Error: An input data type must be specified.' | tee -a "$log_file" >&2
    exit 1
elif [ "$data_type" = "FASTQ" ] || [ "$data_type" = "BAM" ] || [ "$data_type" = "GTF" ]; then
    echo "Input data type: $data_type" | tee -a "$log_file"
else
    echo "Error: The input data type '$data_type' is not supported. Please specify 'FASTQ', 'BAM' or 'GTF'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$organism" ]; then
    echo 'Error: An organism must be specified.' | tee -a "$log_file" >&2
    exit 1
elif [ "$organism" = "HUMAN" ] || [ "$organism" = "CAEEL" ] || [ "$organism" = "DROME" ] || [ "$organism" = "MOUSE" ] || [ "$organism" = "RAT" ] || [ "$organism" = "DANRE" ]; then
    echo "Organism: $organism" | tee -a "$log_file"

    # set the reference protein database according to the organism specified
    if [ "$organism" = "HUMAN" ]; then
        ref_proteome="./data/openprot_uniprotDb_human.txt"
    elif [ "$organism" = "CAEEL" ]; then
        ref_proteome="./data/openprot_uniprotDb_c_elegans.txt"
    elif [ "$organism" = "DROME" ]; then
        ref_proteome="./data/openprot_uniprotDb_drosophila.txt"
    elif [ "$organism" = "MOUSE" ]; then
        ref_proteome="./data/openprot_uniprotDb_mouse.txt"
    elif [ "$organism" = "RAT" ]; then
        ref_proteome="./data/openprot_uniprotDb_rat.txt"
    elif [ "$organism" = "DANRE" ]; then
        ref_proteome="./data/openprot_uniprotDb_zebrafish.txt"
    #elif [ "$organism" = "PANTR" ]; then
    #    ref_proteome="./data/openprot_uniprotDb_chimp.txt"
    #elif [ "$organism" = "BOVIN" ]; then
    #    ref_proteome="./data/openprot_uniprotDb_cow.txt"
    #elif [ "$organism" = "XENTR" ]; then
    #    ref_proteome="./data/openprot_uniprotDb_clawed_frog.txt"
    #elif [ "$organism" = "YEAST" ]; then
    #    ref_proteome="./data/openprot_uniprotDb_yeast.txt"
    fi
else
    echo "Error: The organism '$organism' is not supported. Please specify 'HUMAN', 'CAEEL', 'DROME', 'MOUSE', 'RAT', or 'DANRE'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$reference_gtf" ]; then
    echo 'Error: A reference transcriptome GTF must be specified.' | tee -a "$log_file" >&2
    exit 1
elif ! [ -f "$reference_gtf" ]; then
    echo "Error: The reference transcriptome GTF '$reference_gtf' does not exist but is required." | tee -a "$log_file" >&2
    exit 1
elif ! [ -s "$reference_gtf" ]; then
    echo "Error: The reference transcriptome GTF '$reference_gtf' is empty." | tee -a "$log_file" >&2
    exit 1
elif ! ([[ "$reference_gtf" == *.gtf ]] || [[ "$reference_gtf" == *.gff ]] || [[ "$reference_gtf" == *.gff2 ]] || [[ "$reference_gtf" == *.gff3 ]]); then
    echo "Error: The reference transcriptome GTF '$reference_gtf' must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'." | tee -a "$log_file" >&2
    exit 1
else
    echo "Reference transcriptome GTF: $reference_gtf" | tee -a "$log_file"
fi

if [ -z "$vcf" ]; then
    vcf="None"
    echo 'No VCF provided.' | tee -a "$log_file"
elif ! [ -f "$vcf" ]; then
    echo "Error: The VCF file '$vcf' does not exist." | tee -a "$log_file" >&2
    exit 1
elif ! [ -s "$vcf" ]; then
    echo "Error: The VCF file '$vcf' is empty." | tee -a "$log_file" >&2
    exit 1
elif ! [[ "$vcf" == *.vcf ]]; then
    echo "Error: The VCF file '$vcf' must have a file extension of '.vcf'." | tee -a "$log_file" >&2
    exit 1
else
    echo "VCF file: $vcf" | tee -a "$log_file"
fi

if [ -z "$orf_type" ]; then
    orf_type=all
fi

if [ "$orf_type" = "canonical" ] || [ "$orf_type" = "all" ]; then
    echo "ORF type: $orf_type" | tee -a "$log_file"
else
    echo "Error: The ORF type '$orf_type' is not supported. Please specify either 'canonical' or 'all'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$upstream_orf" ]; then
    upstream_orf=FALSE
fi

if [ "$upstream_orf" = "TRUE" ] || [ "$upstream_orf" = "FALSE" ]; then
    echo "Upstream ORF detection: $upstream_orf" | tee -a "$log_file"
else
    echo "Error: Incorrect boolean value for upstream ORF detection. Please specify either 'TRUE' or 'FALSE'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$downstream_orf" ]; then
    downstream_orf=FALSE
fi

if [ "$downstream_orf" = "TRUE" ] || [ "$downstream_orf" = "FALSE" ]; then
    echo "Downstream ORF detection: $downstream_orf" | tee -a "$log_file"
else
    echo "Error: Incorrect boolean value for downstream ORF detection. Please provide either 'TRUE' or 'FALSE'." | tee -a "$log_file" >&2
    exit 1
fi

if [ -z "$min_orf_length" ]; then
    min_orf_length=30
elif ! is_uint "$min_orf_length"; then
    echo "Error: The minimum ORF length must be a non-negative integer. Entered value: '$min_orf_length'" | tee -a "$log_file" >&2
    exit 1
else
    min_orf_length=$(( 10#"$min_orf_length" ))
    if [ "$min_orf_length" -lt 0 ]; then
        echo "Error: The minimum ORF length must be a non-negative integer. Parsed value: $min_orf_length" | tee -a "$log_file" >&2
        exit 1
    fi
fi

echo "Minimum ORF length (in amino acids): $min_orf_length" | tee -a "$log_file"

if [ -z "$threads" ]; then
    threads=1
elif ! is_uint "$threads"; then
    echo "Error: The number of threads must be a positive integer. Entered value: '$threads'" | tee -a "$log_file" >&2
    exit 1
else
    threads=$(( 10#"$threads" ))
    if [ "$threads" -lt 1 ]; then
        echo "Error: The number of threads must be a positive integer. Parsed value: $threads" | tee -a "$log_file" >&2
        exit 1
    fi
fi

echo "Number of threads: $threads" | tee -a "$log_file"

if [ -z "$min_tx_count" ]; then
    min_tx_count=5
elif ! is_uint "$min_tx_count"; then
    echo "Error: The minimum transcript expression cutoff must be a non-negative integer. Entered value: '$min_tx_count'" | tee -a "$log_file" >&2
    exit 1
else
    min_tx_count=$(( 10#"$min_tx_count" ))
    if [ "$min_tx_count" -lt 0 ]; then
        echo "Error: The minimum transcript expression cutoff must be a non-negative integer. Parsed value: $min_tx_count" | tee -a "$log_file" >&2
        exit 1
    fi
fi

fastq_r1_pe=()
fastq_r2_pe=()
fastq_se=()
bamfiles=()

if [ "$data_type" = "FASTQ" ] || [ "$data_type" = "BAM" ]; then
    if [ -z "$genome_fa" ]; then
        echo "Error: A reference genome FASTA is required for $data_type input files." | tee -a "$log_file" >&2
        exit 1
    elif ! [ -f "$genome_fa" ]; then
        echo "Error: The reference genome FASTA '$genome_fa' does not exist but is required for $data_type input files." | tee -a "$log_file" >&2
        exit 1
    elif ! [ -s "$genome_fa" ]; then
        echo "Error: The reference genome FASTA '$genome_fa' is empty." | tee -a "$log_file" >&2
        exit 1
    elif ! ([[ "$genome_fa" == *.fasta ]] || [[ "$genome_fa" == *.fas ]] || [[ "$genome_fa" == *.fa ]] || [[ "$genome_fa" == *.fna ]] || [[ "$genome_fa" == *.ffn ]] || [[ "$genome_fa" == *.faa ]] || [[ "$genome_fa" == *.mpfa ]] || [[ "$genome_fa" == *.frn ]] || [[ "$genome_fa" == *.fasta.gz ]] || [[ "$genome_fa" == *.fas.gz ]] || [[ "$genome_fa" == *.fa.gz ]] || [[ "$genome_fa" == *.fna.gz ]] || [[ "$genome_fa" == *.ffn.gz ]] || [[ "$genome_fa" == *.faa.gz ]] || [[ "$genome_fa" == *.mpfa.gz ]] || [[ "$genome_fa" == *.frn.gz ]]); then
        echo "Error: The reference genome FASTA '$genome_fa' must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped." | tee -a "$log_file" >&2
        exit 1
    fi

    if [ "$sequencing_platform" = "short-read" ] && [ "$data_type" = "FASTQ" ]; then
        if [ -z "$transcriptome_fa" ]; then
            echo "Error: A reference transcriptome FASTA is required for short-read FASTQ input files." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -f "$transcriptome_fa" ]; then
            echo "Error: The reference transcriptome FASTA '$transcriptome_fa' does not exist but is required for short-read FASTQ input files." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -s "$transcriptome_fa" ]; then
            echo "Error: The reference transcriptome FASTA '$transcriptome_fa' is empty." | tee -a "$log_file" >&2
            exit 1
        elif ! ([[ "$transcriptome_fa" == *.fasta ]] || [[ "$transcriptome_fa" == *.fas ]] || [[ "$transcriptome_fa" == *.fa ]] || [[ "$transcriptome_fa" == *.fna ]] || [[ "$transcriptome_fa" == *.ffn ]] || [[ "$transcriptome_fa" == *.faa ]] || [[ "$transcriptome_fa" == *.mpfa ]] || [[ "$transcriptome_fa" == *.frn ]] || [[ "$transcriptome_fa" == *.fasta.gz ]] || [[ "$transcriptome_fa" == *.fas.gz ]] || [[ "$transcriptome_fa" == *.fa.gz ]] || [[ "$transcriptome_fa" == *.fna.gz ]] || [[ "$transcriptome_fa" == *.ffn.gz ]] || [[ "$transcriptome_fa" == *.faa.gz ]] || [[ "$transcriptome_fa" == *.mpfa.gz ]] || [[ "$transcriptome_fa" == *.frn.gz ]]); then
            echo "Error: The reference transcriptome FASTA '$transcriptome_fa' must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped." | tee -a "$log_file" >&2
            exit 1
        fi
    fi

    if [ -z "$sample_directory" ]; then
        echo "Error: A sample directory must be specified for $data_type input files." | tee -a "$log_file" >&2
        exit 1
    fi

    if [[ ! "$sample_directory" =~ /$ ]]; then
        sample_directory="$sample_directory/"
    fi

    if ! [ -d "$sample_directory" ]; then
        echo "Error: The $data_type sample directory '$sample_directory' does not exist." | tee -a "$log_file" >&2
        exit 1
    fi

    if [ "$data_type" = "FASTQ" ]; then
        # Find all files in the sample directory that have a file extension of '.fastq', '.fq', '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn' (+ '.gz' if gzipped)
        readarray -d '' fastq_file_arr < <(find "$sample_directory" -maxdepth 1 -type f -size +0 -regextype posix-extended -regex ".*\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$" -print0)

        fastq_r1_sample_ids=()
        fastq_r1_files=()
        r1_index=0

        fastq_r2_sample_ids=()
        fastq_r2_files=()
        r2_index=0

        fastq_se_sample_ids=()
        se_index=0

        for file in "${fastq_file_arr[@]}"; do
            filename=$(basename "$file")

            # For short-read FASTQ samples, separate the files into paired-end and single-end
            if [ "$sequencing_platform" = "short-read" ]; then
                sample_id=$(echo "$filename" | sed -E 's/_R1\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
                if [ -n "$sample_id" ] && [ "${#sample_id}" -lt "${#filename}" ]; then
                    fastq_r1_sample_ids[$r1_index]="$sample_id"
                    fastq_r1_files[$r1_index]="$file"
                    r1_index=$(( r1_index + 1 ))
                    continue
                fi

                sample_id=$(echo "$filename" | sed -E 's/_R2\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
                if [ -n "$sample_id" ] && [ "${#sample_id}" -lt "${#filename}" ]; then
                    fastq_r2_sample_ids[$r2_index]="$sample_id"
                    fastq_r2_files[$r2_index]="$file"
                    r2_index=$(( r2_index + 1 ))
                    continue
                fi

                sample_id=$(echo "$filename" | sed -E 's/\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
                if [ -n "$sample_id" ] && [ "${#sample_id}" -lt "${#filename}" ]; then
                    fastq_se_sample_ids[$se_index]="$sample_id"
                    fastq_se[$se_index]="$file"
                    se_index=$(( se_index + 1 ))
                fi

            # For long-read FASTQ samples, do not separate the files into paired-end or single-end
            else
                sample_id=$(echo "$filename" | sed -E 's/\.(fastq|fq|fasta|fas|fa|fna|ffn|faa|mpfa|frn)(\.gz)?$//g')
                if [ -n "$sample_id" ]; then
                    fastq_se_sample_ids[$se_index]="$sample_id"
                    fastq_se[$se_index]="$file"
                    se_index=$(( se_index + 1 ))
                fi
            fi
        done

        # Check for duplicate (single-end) sample IDs
        if [ "$se_index" -gt 1 ]; then
            for index0 in $(seq 0 $(( se_index - 2 ))); do
                sample_id0=${fastq_se_sample_ids[$index0]}
                for index1 in $(seq $(( index0 + 1 )) $(( se_index - 1 ))); do
                    sample_id1=${fastq_se_sample_ids[$index1]}
                    if [ "$sample_id0" = "$sample_id1" ]; then
                        if [ "$sequencing_platform" = "short-read" ]; then
                            echo 'Error: There are duplicate sample IDs for single-end FASTQ files. Please ensure all sample IDs are unique.' | tee -a "$log_file" >&2
                        else
                            echo 'Error: There are duplicate sample IDs for long-read FASTQ files. Please ensure all sample IDs are unique.' | tee -a "$log_file" >&2
                        fi
                        exit 1
                    fi
                done
            done
        fi

        if [ "$sequencing_platform" = "short-read" ]; then
            # Check for duplicate R1 IDs
            for index0 in $(seq 0 $(( r1_index - 2 ))); do
                sample_id0=${fastq_r1_sample_ids[$index0]}
                for index1 in $(seq $(( index0 + 1 )) $(( r1_index - 1 ))); do
                    sample_id1=${fastq_r1_sample_ids[$index1]}
                    if [ "$sample_id0" = "$sample_id1" ]; then
                        echo 'Error: There are duplicate sample IDs for paired-end R1 FASTQ files. Please ensure all sample IDs are unique.' | tee -a "$log_file" >&2
                        exit 1
                    fi
                done
            done

            # Check for duplicate R2 IDs
            for index0 in $(seq 0 $(( r2_index - 2 ))); do
                sample_id0=${fastq_r2_sample_ids[$index0]}
                for index1 in $(seq $(( index0 + 1 )) $(( r2_index - 1 ))); do
                    sample_id1=${fastq_r2_sample_ids[$index1]}
                    if [ "$sample_id0" = "$sample_id1" ]; then
                        echo 'Error: There are duplicate sample IDs for paired-end R2 FASTQ files. Please ensure all sample IDs are unique.' | tee -a "$log_file" >&2
                        exit 1
                    fi
                done
            done

            # Only consider the paired-end FASTQ files that have both ends (R1 and R2)
            # Also make sure all paired-end sample IDs are different from those for single-end
            paired_file_index=0
            for index0 in $(seq 0 $(( r1_index - 1 ))); do
                r1_sample_id=${fastq_r1_sample_ids[$index0]}
                for index1 in $(seq 0 $(( r2_index - 1 ))); do
                    r2_sample_id=${fastq_r2_sample_ids[$index1]}
                    if [ "$r1_sample_id" = "$r2_sample_id" ]; then
                        for se_sample_id in "${fastq_se_sample_ids[@]}"; do
                            if [ "$r1_sample_id" = "$se_sample_id" ]; then
                                echo 'Error: There are matching sample IDs between paired-end and single-end FASTQ files. Please ensure all sample IDs are unique.' | tee -a "$log_file" >&2
                                exit 1
                            fi
                        done

                        fastq_r1_pe[$paired_file_index]=${fastq_r1_files[$index0]}
                        fastq_r2_pe[$paired_file_index]=${fastq_r2_files[$index1]}
                        paired_file_index=$(( paired_file_index + 1 ))
                        break
                    fi
                done
            done

            if [[ ${#fastq_se[@]} -gt 0 ]]; then
                if [ "$r1_index" -gt 0 ] && [ "$r2_index" -eq 0 ]; then
                    echo 'Warning: Paired-end R1 FASTQ samples were found but no R2 FASTQ samples were provided.' | tee -a "$log_file" >&2
                elif [ "$r2_index" -gt 0 ] && [ "$r1_index" -eq 0 ]; then
                    echo 'Warning: Paired-end R2 FASTQ samples were found but no R1 FASTQ samples were provided.' | tee -a "$log_file" >&2
                fi
            elif [[ ${#fastq_r1_pe[@]} -eq 0 ]]; then
                echo "Error: No paired-end or single-end short-read FASTQ samples found in the sample directory '$sample_directory'." | tee -a "$log_file" >&2
                exit 1
            fi
        else
            if [[ ${#fastq_se[@]} -eq 0 ]]; then
                echo "Error: No long-read FASTQ samples found in the sample directory '$sample_directory'." | tee -a "$log_file" >&2
                exit 1
            fi
        fi
    else
        # Find all files in the sample directory that have a file extension of '.bam'
        readarray -d '' bam_file_arr < <(find "$sample_directory" -maxdepth 1 -type f -size +0 -regextype posix-extended -regex ".*\.bam$" -print0)

        bam_file_index=0
        for file in "${bam_file_arr[@]}"; do
            bamfiles[$bam_file_index]="$file"
            bam_file_index=$(( bam_file_index + 1 ))
        done

        if [[ ${#bamfiles[@]} -eq 0 ]]; then
            echo "Error: No BAM samples found in the sample directory '$sample_directory'." | tee -a "$log_file" >&2
            exit 1
        fi
    fi

    echo "$data_type sample directory: $sample_directory" | tee -a "$log_file"

    echo "Reference genome FASTA: $genome_fa" | tee -a "$log_file"

    if [ "$sequencing_platform" = "short-read" ] && [ "$data_type" = "FASTQ" ]; then
        echo "Reference transcriptome FASTA: $transcriptome_fa" | tee -a "$log_file"
    fi
else
    if [ "$sequencing_platform" = "long-read" ]; then
        if [ -z "$custom_gtf" ]; then
            echo "Error: A custom transcript annotation GTF is required for long-read GTF input." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -f "$custom_gtf" ]; then
            echo "Error: The custom transcript annotation GTF '$custom_gtf' does not exist but is required for long-read GTF input." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -s "$custom_gtf" ]; then
            echo "Error: The custom transcript annotation GTF '$custom_gtf' is empty." | tee -a "$log_file" >&2
            exit 1
        elif ! ([[ "$custom_gtf" == *.gtf ]] || [[ "$custom_gtf" == *.gff ]] || [[ "$custom_gtf" == *.gff2 ]] || [[ "$custom_gtf" == *.gff3 ]]); then
            echo "Error: The custom transcript annotation GTF '$custom_gtf' must have a file extension of '.gtf', '.gff', '.gff2' or '.gff3'." | tee -a "$log_file" >&2
            exit 1
        fi

        echo "Custom transcript annotation GTF: $custom_gtf" | tee -a "$log_file"
    fi

    if [ "$vcf" != "None" ]; then
        if [ -z "$genome_fa" ]; then
            echo "Error: A reference genome FASTA is required for generating a variant-aware proteome database." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -f "$genome_fa" ]; then
            echo "Error: The reference genome FASTA '$genome_fa' does not exist but is required for generating a variant-aware proteome database." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -s "$genome_fa" ]; then
            echo "Error: The reference genome FASTA '$genome_fa' is empty." | tee -a "$log_file" >&2
            exit 1
        elif ! ([[ "$genome_fa" == *.fasta ]] || [[ "$genome_fa" == *.fas ]] || [[ "$genome_fa" == *.fa ]] || [[ "$genome_fa" == *.fna ]] || [[ "$genome_fa" == *.ffn ]] || [[ "$genome_fa" == *.faa ]] || [[ "$genome_fa" == *.mpfa ]] || [[ "$genome_fa" == *.frn ]] || [[ "$genome_fa" == *.fasta.gz ]] || [[ "$genome_fa" == *.fas.gz ]] || [[ "$genome_fa" == *.fa.gz ]] || [[ "$genome_fa" == *.fna.gz ]] || [[ "$genome_fa" == *.ffn.gz ]] || [[ "$genome_fa" == *.faa.gz ]] || [[ "$genome_fa" == *.mpfa.gz ]] || [[ "$genome_fa" == *.frn.gz ]]); then
            echo "Error: The reference genome FASTA '$genome_fa' must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn', or any of these file extensions with '.gz' appended if the file is gzipped." | tee -a "$log_file" >&2
            exit 1
        fi

        echo "Reference genome FASTA: $genome_fa" | tee -a "$log_file"
    fi

    if [ "$sequencing_platform" = "long-read" ]; then
        if [ -z "$tx_counts_file" ]; then
            tx_counts_file="None"
        elif ! [ -f "$tx_counts_file" ]; then
            echo "Error: The transcript counts file '$tx_counts_file' does not exist." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -s "$tx_counts_file" ]; then
            echo "Error: The transcript counts file '$tx_counts_file' is empty." | tee -a "$log_file" >&2
            exit 1
        elif ! ([[ "$tx_counts_file" == *.txt ]] || [[ "$tx_counts_file" == *.csv ]] || [[ "$custom_gtf" == *.tsv ]]); then
            echo "Error: The transcript counts file '$tx_counts_file' must have a file extension of '.txt', '.csv' or '.tsv'." | tee -a "$log_file" >&2
            exit 1
        fi
    else
        if [ -z "$tx_counts_file" ]; then
            echo "Error: A transcript counts file is required for short-read GTF input." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -f "$tx_counts_file" ]; then
            echo "Error: The transcript counts file '$tx_counts_file' does not exist but is required for short-read GTF input." | tee -a "$log_file" >&2
            exit 1
        elif ! [ -s "$tx_counts_file" ]; then
            echo "Error: The transcript counts file '$tx_counts_file' is empty." | tee -a "$log_file" >&2
            exit 1
        elif ! ([[ "$tx_counts_file" == *.txt ]] || [[ "$tx_counts_file" == *.csv ]] || [[ "$custom_gtf" == *.tsv ]]); then
            echo "Error: The transcript counts file '$tx_counts_file' must have a file extension of '.txt', '.csv' or '.tsv'." | tee -a "$log_file" >&2
            exit 1
        fi
    fi
fi

# Capture the time the run started
start=$(date)
start_time=$(date +%s)
echo "Run started at: $start" | tee -a "$log_file"

if [ "$data_type" = "GTF" ]; then
    custom_gtf_file="$custom_gtf"
    if [ "$sequencing_platform" = "long-read" ]; then
        custom_gtf_file="$reference_gtf"
    fi
elif [ "$sequencing_platform" = "short-read" ]; then
    process_short_read_samples fastq_r1_pe[@] fastq_r2_pe[@] fastq_se[@] bamfiles[@] "$genome_fa" "$transcriptome_fa" "$reference_gtf" "$output_directory" "$threads"

    custom_gtf_file="$reference_gtf"
    tx_counts_file="${output_directory}salmon_transcript_counts.txt"
else
    process_long_read_samples fastq_se[@] "$genome_fa" "$reference_gtf" "$sample_directory" "$output_directory" "$threads" "$organism"

    custom_gtf_file="${output_directory}bambu_transcript_annotations.gtf"
    tx_counts_file="${output_directory}bambu_transcript_counts.txt"
fi

ref_genome_file="$genome_fa"

# If the reference genome is gzipped, decompress it first
if [[ "$genome_fa" == *.gz ]]; then
    temp_genome_fa="${output_directory}decompressed_genome.fa"
    if ! [ -f "$temp_genome_fa" ]; then
        gzip -c -d "$ref_genome_file" > "$temp_genome_fa"
    fi
    ref_genome_file="$temp_genome_fa"
fi

# Generate the proteome database
generate_database "$custom_gtf_file" "$output_directory" "$reference_gtf" "$organism" "$min_orf_length" "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$ref_genome_file" "$tx_counts_file" "$min_tx_count" "$data_type" "$threads"

# Delete the temporary decompressed reference genome
if [[ "$genome_fa" == *.gz ]]; then
    rm "$ref_genome_file"
fi

# Capture the time the run ended
end=$(date)
end_time=$(date +%s)
echo "Run ended at: $end" | tee -a "$log_file"

# Convert the runtime to days, hours, minutes and seconds
duration=$(( end_time - start_time ))
days=$(( duration / 86400 ))
hours=$(( (duration % 86400) / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$(( duration % 60 ))
echo "Runtime: $days days, $hours hours, $minutes minutes, $seconds seconds" | tee -a "$log_file"