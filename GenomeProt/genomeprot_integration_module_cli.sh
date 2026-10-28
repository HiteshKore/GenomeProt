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

function reformat_proteomics_results() {
    input_dir=$1
    search_tool=$2
    reformatted_file=$3

    echo "Reformatting $search_tool proteomics results in '$input_dir'" | tee -a "$log_file"
    "$conda_path" run -n "$genomeprot_env" Rscript ./bin/integration_module/reformat_peptide_results.R -d "$input_dir" -s "$search_tool" -o "$reformatted_file"
    echo "Reformatted proteomics results written to '$reformatted_file'." | tee -a "$log_file"
}

function integrate_proteogenomics() {
    proteomics=$1
    metadata=$2
    gtf=$3
    output_dir=$4

    echo "Reformatted proteomics results file: $proteomics"
    echo "Proteome database metadata file: $metadata"
    echo "Proteome database transcripts GTF file: $gtf"
    echo "Output directory: $output_dir"

    echo 'Performing proteogenomics integration' | tee -a "$log_file"
    "$conda_path" run -n "$genomeprot_env" Rscript ./bin/integration_module/map_peptides_generate_outputs.R -p "$proteomics" -m "$metadata" -g "$gtf" -s "$output_dir"

    echo 'Generating proteogenomics summary report' | tee -a "$log_file"
    "$conda_path" run -n "$genomeprot_env" Rscript ./bin/integration_module/generate_integration_summary_report.R -i ./bin/integration_module/integration_summary_report.Rmd -o "$output_dir"

    echo 'Proteogenomics integration completed.' | tee -a "$log_file"
}

function usage() {
    cat << EOF
Program: GenomeProt proteogenomics integration module (command-line)
Version: 0.0.1v

 ▗▄▄▖▗▄▄▄▖▗▖  ▗▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▖ ▗▄▄▖  ▗▄▖▗▄▄▄▖
▐▌   ▐▌   ▐▛▚▖▐▌▐▌ ▐▌▐▛▚▞▜▌▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌ █
▐▌▝▜▌▐▛▀▀▘▐▌ ▝▜▌▐▌ ▐▌▐▌  ▐▌▐▛▀▀▘▐▛▀▘ ▐▛▀▚▖▐▌ ▐▌ █
▝▚▄▞▘▐▙▄▄▖▐▌  ▐▌▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▌   ▐▌ ▐▌▝▚▄▞▘ █

Usage: genomeprot_integration_module_cli.sh [-h] [-d <proteomics results dir>] [-t <proteomics search tool>] [-r <path to write reformatted proteomics results>] [-p <reformatted proteomics results file>] [-m <proteome database metadata file>] [-g <proteome database transcripts GTF file>] [-o <output dir>]

Integrate transcriptomics and proteomics results.

GENERAL
  -h, --help                            Display this help and exit

1) PROTEOMICS RESULTS REFORMATTING
  -d, --proteomics-results-dir <dir>    Proteomics results directory; must contain .txt, .csv or .tsv proteomics results files
  -t, --search-tool <tool>              Search tool used to generate the proteomic data: Spectronaut | FragPipe | FragPipe_quant | DIA-NN | MaxQuant
  -r, --reformatted-output <path>       Path to write the reformatted proteomics results file to

2) PROTEOGENOMICS INTEGRATION
  -p, --proteomics <file>               Reformatted proteomics results file
  -m, --metadata <file>                 Proteome database metadata file (e.g. proteome_database_metadata.txt)
  -g, --gtf <file>                      Proteome database transcripts GTF file (e.g. proteome_database_transcripts.gtf)
  -o, --output-dir <dir>                Output directory

=== Modes ===

This script can 1) reformat proteomics results and 2) perform proteogenomics integration. Users can choose to run both modes or only one mode.

To run both modes, specify the -d, -t, -r, -m, -g and -o arguments.
To run one mode only, specify all arguments for that mode.

=== Proteomics search tools ===

Spectronaut: Spectronaut output data.
FragPipe: peptide.tsv from FragPipe.
FragPipe_quant: report.pr_matrix.tsv from FragPipe.
DIA-NN: report.pr_matrix.tsv from DIA-NN.
MaxQuant: peptides.txt from MaxQuant.

=== Expected file formats ===

All proteomics results files in the proteomics results directory (-d / --proteomics-results-dir) must have a file extension of '.txt', '.csv' or '.tsv', and must be named as '<dataset ID><file extension>', where all dataset IDs are unique and non-empty.
The reformatted proteomics results file (-p / --proteomics) must have a file extension of '.txt', '.csv' or '.tsv'.
The proteome database metadata file (-m / --metadata) must have a file extension of '.txt'.
The proteome database transcripts GTF file (-g / --gtf) must have a file extension of '.gtf'.

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
OPTS=$(getopt -o d:t:r:p:m:g:o:h \
              -l proteomics-results-dir:,search-tool:,reformatted-output:,proteomics:,metadata:,gtf:,output-dir:,help \
              -n "genomeprot_integration_module_cli.sh" -- "$@")

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
        -h | --help )                   usage; exit 0 ;;
        -d | --proteomics-results-dir ) proteomics_directory="$2"; shift 2 ;;
        -t | --search-tool )            search_tool="$2"; shift 2 ;;
        -r | --reformatted-output )     reformatted_file="$2"; shift 2 ;;
        -p | --proteomics )             proteomics="$2"; shift 2 ;;
        -m | --metadata )               metadata="$2"; shift 2 ;;
        -g | --gtf )                    gtf="$2"; shift 2 ;;
        -o | --output-dir )             output_directory="$2"; shift 2 ;;
        -- ) shift; break ;;        # End of options
        * ) echo "Invalid option: $1"; usage; exit 1 ;;
    esac
done

# Check whether the input arguments are correct

if [ -n "$proteomics_directory" ]; then
    if [[ ! "$proteomics_directory" =~ /$ ]]; then
        proteomics_directory="$proteomics_directory/"
    fi

    if ! [ -d "$proteomics_directory" ]; then
        echo "Error: The proteomics results directory '$proteomics_directory' does not exist." >&2
        exit 1
    fi

    # Find all files in the proteomics results directory that have a file extension of '.txt', '.csv' or '.tsv'
    readarray -d '' proteomics_file_arr < <(find "$proteomics_directory" -maxdepth 1 -type f -size +0 -regextype posix-extended -regex ".*\.(txt|csv|tsv)$" -print0)

    dataset_ids=()
    dataset_id_index=0

    # Get the dataset ID for each file found
    for file in "${proteomics_file_arr[@]}"; do
        filename=$(basename "$file")
        dataset_id=$(echo "$filename" | sed -E 's/\.(txt|csv|tsv)$//g')
        if [ -z "$dataset_id" ]; then
            echo "Error: There are empty dataset IDs in the proteomics results directory '$proteomics_directory'. Please ensure all dataset IDs consist of at least one character." >&2
            exit 1
        fi
        dataset_ids[$dataset_id_index]="$dataset_id"
        dataset_id_index=$(( dataset_id_index + 1 ))
    done

    if [ "$dataset_id_index" -gt 1 ]; then
        for index0 in $(seq 0 $(( dataset_id_index - 2 ))); do
            dataset_id0=${dataset_ids[$index0]}
            for index1 in $(seq $(( index0 + 1 )) $(( dataset_id_index - 1 ))); do
                dataset_id1=${dataset_ids[$index1]}
                if [ "$dataset_id0" = "$dataset_id1" ]; then
                    echo "Error: There are duplicate dataset IDs in the proteomics results directory '$proteomics_directory'. Please ensure all dataset IDs are unique." >&2
                    exit 1
                fi
            done
        done
    elif [ "$dataset_id_index" -eq 0 ]; then
        echo 'Error: Please provide non-empty proteomics results files in the proteomics results directory specified.' >&2
        exit 1
    fi
fi

if [ -n "$search_tool" ] && [ "$search_tool" != "Spectronaut" ] && [ "$search_tool" != "FragPipe" ] && [ "$search_tool" != "FragPipe_quant" ] && [ "$search_tool" != "DIA-NN" ] && [ "$search_tool" != "MaxQuant" ]; then
    echo "Error: The proteomics search tool '$search_tool' is not supported. Please specify 'Spectronaut', 'FragPipe', 'FragPipe_quant', 'DIA-NN' or 'MaxQuant'." >&2
    exit 1
fi

# Create the directory for storing the reformatted peptide results file if it does not exist
if [ -n "$reformatted_file" ]; then
    reformatted_file_dir=$(dirname "$reformatted_file")
    if ! [ -d "$reformatted_file_dir" ]; then
        mkdir -p "$reformatted_file_dir"
    fi
fi

if [ -z "$reformatted_file" ] && [ -n "$proteomics" ]; then
    if ! [ -f "$proteomics" ]; then
        echo "Error: The reformatted proteomics results file '$proteomics' does not exist." >&2
        exit 1
    elif ! [ -s "$proteomics" ]; then
        echo "Error: The reformatted proteomics results file '$proteomics' is empty." >&2
        exit 1
    elif ! ([[ "$proteomics" == *.txt ]] || [[ "$proteomics" == *.csv ]] || [[ "$proteomics" == *.tsv ]]); then
        echo "Error: The reformatted proteomics results file '$proteomics' must have a file extension of '.txt', '.csv' or '.tsv'." >&2
        exit 1
    fi
fi

if [ -n "$metadata" ]; then
    if ! [ -f "$metadata" ]; then
        echo "Error: The proteome database metadata file '$metadata' does not exist." >&2
        exit 1
    elif ! [ -s "$metadata" ]; then
        echo "Error: The proteome database metadata file '$metadata' is empty." >&2
        exit 1
    elif ! [[ "$metadata" == *.txt ]]; then
        echo "Error: The proteome database metadata file '$metadata' must have a file extension of '.txt'." >&2
        exit 1
    fi
fi

if [ -n "$gtf" ]; then
    if ! [ -f "$gtf" ]; then
        echo "Error: The proteome database transcripts GTF file '$gtf' does not exist." >&2
        exit 1
    elif ! [ -s "$gtf" ]; then
        echo "Error: The proteome database transcripts GTF file '$gtf' is empty." >&2
        exit 1
    elif ! [[ "$gtf" == *.gtf ]]; then
        echo "Error: The proteome database transcripts GTF file '$gtf' must have a file extension of '.gtf'." >&2
        exit 1
    fi
fi

# Create the output directory if it does not exist
if [ -n "$output_directory" ]; then
    if [[ ! "$output_directory" =~ /$ ]]; then
        output_directory="$output_directory/"
    fi

    if ! [ -d "$output_directory" ]; then
        mkdir -p "$output_directory"
    fi

    log_file="${output_directory}integration_module_logfile.txt"
else
    log_file=integration_module_logfile.txt
fi

if [ -n "$proteomics_directory" ] && [ -n "$search_tool" ] && [ -n "$reformatted_file" ] && [ -n "$metadata" ] && [ -n "$gtf" ] && [ -n "$output_directory" ]; then
    run_both_modes=true
elif [ -n "$proteomics_directory" ] && [ -n "$search_tool" ] && [ -n "$reformatted_file" ]; then
    run_reformat_mode_only=true
elif [ -n "$proteomics" ] && [ -n "$metadata" ] && [ -n "$gtf" ] && [ -n "$output_directory" ]; then
    run_reformat_mode_only=false
else
    echo 'Error: Please specify sufficient arguments to run both modes (proteomics results reformatting and proteogenomics integration) or only one mode.' >&2
    exit 1
fi

echo "Log file: $log_file" | tee -a "$log_file"

# Capture the time the run started
start=$(date)
start_time=$(date +%s)
echo "Run started at: $start" | tee -a "$log_file"

if [ "$run_both_modes" = "true" ]; then
    reformat_proteomics_results "$proteomics_directory" "$search_tool" "$reformatted_file"
    integrate_proteogenomics "$reformatted_file" "$metadata" "$gtf" "$output_directory"
elif [ "$run_reformat_mode_only" = "true" ]; then
    reformat_proteomics_results "$proteomics_directory" "$search_tool" "$reformatted_file"
else
    integrate_proteogenomics "$proteomics" "$metadata" "$gtf" "$output_directory"
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