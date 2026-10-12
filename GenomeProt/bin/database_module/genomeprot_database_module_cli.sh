#!/bin/bash

# Get the path to the conda binary
conda_path=$(/usr/bin/which conda)
if [ -z "$conda_path" ]; then
    echo 'Failed to find the conda binary. Check if conda is installed or if it is within a directory specified in the $PATH variable.' >&2
    exit 1
fi

# Name of the Conda environment
genomeprot_env="GenomeProt_env"

# Check if the Conda environment exists
if ! ($conda_path env list | grep -q "^$genomeprot_env\s"); then
    echo "Conda environment '$genomeprot_env' does not exist."
    exit 1
fi

function short_read_commands() {
    fq_se=("${!1}")
    bm_fn=("${!2}")
    fastq_id_pe=("${!3}")
    genome_fa=$4
    transcriptome_db=$5
    ref_gtf=$6
    sample_dir=$7
    output_dir=$8
    organism=$9
    orf_length=${10}
    vcf_fn=${11}
    orf_input_type=${12}
    threads=${13}
    five_utr_orf=${14}
    three_utr_orf=${15}
    ref_proteome=${16}

    
    echo ${bm_fn[@]}
    echo "Generating salmon index"
    if [[ ${#fastq_id_pe[@]} -gt 0 || ${#fq_se[@]} -gt 0 ]]; then
          # FASTQ samples found, so build the salmon index

          # Check if both the reference genome and reference transcriptome are gzipped
          if [[ "$transcriptome_db" == *.gz && "$genome_fa" == *.gz ]]; then
              echo "Processing gzipped files" | tee -a "$log_file"
      
              # Generate decoy file (sequence names without '>')
              gunzip -c "$genome_fa" | grep '^>' | cut -d ' ' -f 1 | sed 's/>//' > "${output_dir}decoys.txt"
      
              # Concatenate transcriptome and genome to create gentrome
              cat "$transcriptome_db" "$genome_fa" > "${output_dir}gentrome.fa.gz"
      
              # Run Salmon index
              "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                  -t "${output_dir}gentrome.fa.gz" -d "${output_dir}decoys.txt" \
                  -p "$threads" -i "${output_dir}salmon_index"

          elif [[ "$transcriptome_db" != *.gz && "$genome_fa" != *.gz ]]; then
              echo "Processing uncompressed files" | tee -a "$log_file"
      
              # Generate decoy file (sequence names without '>')
              grep '^>' "$genome_fa" | cut -d ' ' -f 1 | sed 's/>//' > "${output_dir}decoys.txt"
      
              # Concatenate transcriptome and genome to create gentrome
              cat "$transcriptome_db" "$genome_fa" > "${output_dir}gentrome.fa"
      
              # Run Salmon index
              "$conda_path" run -n "$genomeprot_env" salmon index --gencode \
                  -t "${output_dir}gentrome.fa" -d "${output_dir}decoys.txt" \
                  -p "$threads" -i "${output_dir}salmon_index"

          else
              echo "Error: the transcriptome and genome files must both be gzipped or both uncompressed." | tee -a "$log_file" >&2
              exit 1
          fi

    else
        echo "No FASTQ samples found, skipping salmon index generation." | tee -a "$log_file"
    fi

    #processing files
    if [[ ${#fastq_id_pe[@]} -gt 0 ]]; then
        # Process paired-end samples
        for id in "${fastq_id_pe[@]}"; do
            echo "Quantifying paired-end reads from $id" | tee -a "$log_file"
    
            r1_file_path=$(find "$sample_dir" -type f \( -name "${id}_R1.fastq" -o -name "${id}_R1.fastq.gz" \))
            r2_file_path=$(find "$sample_dir" -type f \( -name "${id}_R2.fastq" -o -name "${id}_R2.fastq.gz" \))
    
            "$conda_path" run -n "$genomeprot_env" salmon quant \
                -i "${output_dir}salmon_index" -p "$threads" -l A \
                -1 "$r1_file_path" -2 "$r2_file_path" \
                -o "${output_dir}${id}"
        done
    
        echo "Paired-end samples processed." | tee -a "$log_file"

    elif [[ ${#fq_se[@]} -gt 0 ]]; then
        # Process single-end samples
        for file in "${fq_se[@]}"; do
            id="${file%.gz}"
            id="${id%.fastq}"
            echo "Quantifying single-end reads from $id" | tee -a "$log_file"
    
            "$conda_path" run -n "$genomeprot_env" salmon quant \
                -i "${output_dir}salmon_index" -p "$threads" -l A \
                -r "${sample_dir}${file}" \
                -o "${output_dir}${id}"
        done
    
        echo "Single-end samples processed." | tee -a "$log_file"

    elif [[ ${#bm_fn[@]} -gt 0 ]]; then
        # Process BAM samples
        for file in "${bm_fn[@]}"; do
            id="${file%.bam}"
            echo "Quantifying reads from $id" | tee -a "$log_file"
    
            "$conda_path" run -n "$genomeprot_env" salmon quant \
                -t "$transcriptome_db" -p "$threads" -l A \
                -a "${sample_dir}${file}" \
                -o "${output_dir}${id}"
        done
    
        echo "BAM samples processed." | tee -a "$log_file"
    
    else
        echo "Error: no paired-end, single-end or BAM samples found in $sample_dir" | tee -a "$log_file" >&2
        exit 1
    fi

    echo "Generating transcript counts from Salmon quantification files" | tee -a "$log_file"
    
    $conda_path run -n $genomeprot_env Rscript ./bin/database_module/matrix_compilation_salmon.R -s "$output_dir" -g "$ref_gtf" -o $output_dir"salmon_transcript_counts.txt"
}

function generate_database() {
    custom_gtf=$1
    output_dir=$2
    ref_gtf=$3
    organism=$4
    orf_length=$5
    vcf_fn=$6
    orf_input_type=$7
    five_utr_orf=$8
    three_utr_orf=$9
    ref_proteome=${10}
    ref_genome=${11}
    tx_count_file=${12}
    min_tx_count=${13}
    input_type=${14}

    if [[ "$vcf_fn" == "None" ]]; then      # if no VCF file was provided
        echo "Generating the proteome database" | tee -a "$log_file"
        echo $tx_count_file
        echo $min_tx_count

        if [[ "$tx_count_file" == "None" ]]; then
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/generate_proteome.R -g $custom_gtf -r $ref_gtf -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -s $output_dir
        else
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/generate_proteome.R -g $custom_gtf -r $ref_gtf -c $tx_count_file -m $min_tx_count -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -s $output_dir
        fi

        echo "Annotating the proteome database" | tee -a "$log_file"

        if [[ "$input_type" == "GTF" ]]; then
            $conda_path run -n $genomeprot_env python ./bin/database_module/annotate_proteome.py $ref_gtf $ref_proteome $output_dir"ORFome_aa.txt" $custom_gtf $output_dir $orf_input_type $orf_length None $organism $threads 2000
        else
            $conda_path run -n $genomeprot_env python ./bin/database_module/annotate_proteome.py $ref_gtf $ref_proteome $output_dir"ORFome_aa.txt" $output_dir"proteome_database_transcripts.gtf" $output_dir $orf_input_type $orf_length None $organism $threads 2000
        fi
    else                                    # if a VCF file was provided
        echo "Generating the proteome database" | tee -a "$log_file"
        echo $tx_count_file
        echo $min_tx_count

        if [[ "$tx_count_file" == "None" ]]; then
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/generate_proteome.R -G $ref_genome -g $custom_gtf -r $ref_gtf -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -v $vcf_fn -s $output_dir
        else
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/generate_proteome.R -G $ref_genome -g $custom_gtf -r $ref_gtf -c $tx_count_file -m $min_tx_count -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -v $vcf_fn -s $output_dir
        fi

        echo "Annotating the proteome database" | tee -a "$log_file"

        mutant_db=$output_dir"Mutant_ORFome_aa.txt"
        $conda_path run -n $genomeprot_env python ./bin/database_module/annotate_proteome.py $ref_gtf $ref_proteome $output_dir"ORFome_aa.txt" $output_dir"proteome_database_transcripts.gtf" $output_dir $orf_input_type $orf_length $mutant_db $organism $threads 2000
    fi
}

function fastq_bam_input() {
    fq_se=("${!1}")    # Use indirect referencing to dereference the array
    bm_fn=("${!2}")
    fastq_id_pe=("${!3}")
    platform=$4
    genome_fa=$5
    sample_dir=$6
    output_dir=$7
    threads=$8
    ref_gtf=$9
    organism=${10}
    orf_length=${11}
    vcf_fn=${12}
    orf_input_type=${13}
    five_utr_orf=${14}
    three_utr_orf=${15}
    ref_proteome=${16}
    transcript_expr_cutoff=${17}
    input_type=${18}

    index_fn=$(basename "$genome_fa")".mmi"

    # currently supports long-read data only
    if [[ $platform == "long-read" ]]; then
        # FASTQ files
        if [[ ${#fq_se[@]} -gt 0 ]]; then
            echo "Generating minimap2 index" | tee -a "$log_file"
            $conda_path run -n $genomeprot_env bash -c "minimap2 -ax splice:hq -d  $output_dir$index_fn  $genome_fa"
            echo "Indexing complete" | tee -a "$log_file"

            echo "Genome alignment started" | tee -a "$log_file"
            for file in "${fq_se[@]}"; do
                file_prefix="${file%%.fastq*}"
                echo "Processing:\t"$file_prefix | tee -a "$log_file"
                $conda_path run -n $genomeprot_env bash -c "minimap2 -t $threads -ax splice:hq --sam-hit-only --secondary=no $output_dir$index_fn $sample_dir$file | samtools view -bh | samtools sort -@ $threads -o $output_dir$file_prefix.bam"
            done
            echo "Genome alignment completed" | tee -a "$log_file"

            echo "Performing transcript identification and quantification with bambu" | tee -a "$log_file"
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/run_bambu.R -b $output_dir -g $ref_gtf -t $threads -s $organism -o $output_dir 2>&1 | tee -a "$log_file"
        fi

        # BAM files
        if [[ ${#bm_fn[@]} -gt 0 ]]; then
            echo "Performing transcript identification and quantification with bambu" | tee -a "$log_file"
            $conda_path run -n $genomeprot_env Rscript ./bin/database_module/run_bambu.R -b $sample_dir -g $ref_gtf -t $threads -s $organism -o $output_dir 2>&1 | tee -a "$log_file"
        fi

        transcript_count_fn=$output_directory"bambu_transcript_counts.txt"
        custom_gtf=$output_directory"bambu_transcript_annotations.gtf"

        #protome database generation
        generate_database "$custom_gtf" "$output_dir" "$ref_gtf" "$organism" "$orf_length" "$vcf_fn" "$orf_input_type" "$five_utr_orf" "$three_utr_orf" "$ref_proteome" "$genome_fa" "$transcript_count_fn" "$transcript_expr_cutoff" "$input_type"
    fi
}

function usage() {
    cat << EOF
Program: GenomeProt Database Module Commandline
Version: 0.0.1v

 ▗▄▄▖▗▄▄▄▖▗▖  ▗▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▖ ▗▄▄▖  ▗▄▖▗▄▄▄▖
▐▌   ▐▌   ▐▛▚▖▐▌▐▌ ▐▌▐▛▚▞▜▌▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌ █
▐▌▝▜▌▐▛▀▀▘▐▌ ▝▜▌▐▌ ▐▌▐▌  ▐▌▐▛▀▀▘▐▛▀▘ ▐▛▀▚▖▐▌ ▐▌ █
▝▚▄▞▘▐▙▄▄▖▐▌  ▐▌▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▌   ▐▌ ▐▌▝▚▄▞▘ █

Usage: genomeprot_db_generation.sh [-h] [-s <sequencing type>] [-o <organism>] [-l <ORF length (AA)>] [-t <data type>] [-g <reference genome fasta>] [-p <threads>] [-a <reference GTF>] [-G <custom GTF>] [-v <VCF file>] [-T <ORF type>] [-r <transcriptome database>] [-U <upstream ORFs>] [-D <downstream ORFs>] [-m <minimum transcript count>] [-C <transcript count file>] [-d <sample directory>] [-O <output directory>]

Generate a proteome database from long-read or short-read RNA-seq data.
Supported inputs: FASTQ, BAM or GTF.

GENERAL
  -h, --help                      Display this help and exit
  -p, --threads <int>             Number of threads (default: 1)
  -O, --output-dir <dir>          Output directory

INPUT
  -s, --sequencing-type <type>     Sequencing type: long-read | short-read
  -t, --data-type <type>           Input data type: FASTQ | BAM | GTF
  -d, --sample-dir <dir>           Sample directory [required for FASTQ and BAM input]
  -o, --organism <name>            Organism: HUMAN | CAEEL | DROME | MOUSE | RAT | DANRE

REFERENCES
  -g, --genome <file>             Reference genome FASTA [required for FASTQ and BAM input]
  -a, --reference-gtf <file>      Reference GTF (GENCODE annotations only)
  -G, --custom-gtf <file>         Custom GTF [required for BAM and GTF input]
  -r, --transcriptome-db <file>   Transcriptome database [required for short-read RNA-Seq FASTQ and BAM input]
  -v, --vcf <file>                VCF file with variants from all samples [optional; use for variant-aware proteome DB]

ORF OPTIONS
  -T, --orf-type <type>           ORF type: canonical | all
  -l, --orf-length <int>          Minimum ORF length in amino acids (default: 30)
  -U, --upstream-orf <bool>       Include upstream (5'UTR) ORFs: TRUE | FALSE
  -D, --downstream-orf <bool>     Include downstream (3'UTR) ORFs: TRUE | FALSE

EXPRESSION FILTERING
  -m, --min-tx-count <int>        Minimum transcript expression cutoff (default: 5) [optional]
  -C, --tx-count-file <file>      Transcript count file [optional]

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
              -l sequencing-type:,organism:,orf-length:,data-type:,genome:,threads:,reference-gtf:,custom-gtf:,vcf:,orf-type:,transcriptome-db:,upstream-orf:,downstream-orf:,min-tx-count:,tx-count-file:,sample-dir:,output-dir:,help \
              -n "genomeprot_db_generation.sh" -- "$@")

if [ $? != 0 ]; then
    echo "Failed to parse options." >&2
    usage
    exit 1
fi

# Re-evaluate the set of positional parameters
eval set -- "$OPTS"

# Process options
while true; do
    case "$1" in
        -s | --sequencing-type )    sequencing_platform="$2"; shift 2 ;;
        -o | --organism )           organism="$2"; shift 2 ;;
        -l | --orf-length )         orf_length="$2"; shift 2 ;;
        -t | --data-type )          data_type="$2"; shift 2 ;;
        -g | --genome )             genome_fa="$2"; shift 2 ;;
        -p | --threads )            threads="$2"; shift 2 ;;
        -a | --reference-gtf )      reference_gtf="$2"; shift 2 ;;
        -G | --custom-gtf )         custom_gtf="$2"; shift 2 ;;
        -v | --vcf )                vcf="$2"; shift 2 ;;
        -T | --orf-type )           orf_type="$2"; shift 2 ;;
        -r | --transcriptome-db )   transcriptome_db="$2"; shift 2 ;;
        -U | --upstream-orf )       upstream_orf="$2"; shift 2 ;;
        -D | --downstream-orf )     downstream_orf="$2"; shift 2 ;;
        -m | --min-tx-count )       min_tx_count="$2"; shift 2 ;;
        -C | --tx-count-file )      tx_count_file="$2"; shift 2 ;;
        -d | --sample-dir )         sample_directory="$2"; shift 2 ;;
        -O | --output-dir )         output_directory="$2"; shift 2 ;;
        -h | --help )               usage; exit 0 ;;  # Call usage function and exit
        -- ) shift; break ;;        # End of options
        * ) echo "Invalid option: $1"; usage; exit 1 ;;
    esac
done

# Check whether the input arguments are correct

if [[ ! "$output_directory" =~ /$ ]]; then
    output_directory="$output_directory/"
fi

if [[ ! "$sample_directory" =~ /$ ]]; then
    sample_directory="$sample_directory/"
fi

# log file
log_file=$output_directory"logfile.txt"

echo $log_file

if [ "$sequencing_platform" = "long-read" ] || [ "$sequencing_platform" = "short-read" ]; then
    echo "Sequencing type: $sequencing_platform" | tee -a "$log_file"
else
    echo "The sequencing type '$sequencing_platform' is not supported. Please specify either 'long-read' or 'short-read'." | tee -a "$log_file"
    exit 1
fi

if [ "$organism" = "HUMAN" ] || [ "$organism" = "CAEEL" ] || [ "$organism" = "DROME" ] || [ "$organism" = "MOUSE" ] || [ "$organism" = "RAT" ] || [ "$organism" = "DANRE" ]; then
    echo "Organism: $organism" | tee -a "$log_file"

    # set the reference protein database according to the organism specified
    if [[ "$organism" == "HUMAN" ]]; then
        ref_proteome="./data/openprot_uniprotDb_human.txt"
    elif [[ "$organism" == "CAEEL" ]]; then
        ref_proteome="./data/openprot_uniprotDb_c_elegans.txt"
    elif [[ "$organism" == "DROME" ]]; then
        ref_proteome="./data/openprot_uniprotDb_drosophila.txt"
    elif [[ "$organism" == "MOUSE" ]]; then
        ref_proteome="./data/openprot_uniprotDb_mouse.txt"
    elif [[ "$organism" == "RAT" ]]; then
        ref_proteome="./data/openprot_uniprotDb_rat.txt"
    elif [[ "$organism" == "DANRE" ]]; then
        ref_proteome="./data/openprot_uniprotDb_zebrafish.txt"
    fi
else
    echo "The organism '$organism' is not supported. Please specify 'HUMAN', 'CAEEL', 'DROME', 'MOUSE', 'RAT', or 'DANRE'." | tee -a "$log_file"
    exit 1
fi

if [ -n "$vcf" ]; then
    echo "VCF file: $vcf" | tee -a "$log_file"
else
    vcf="None"
fi

if [ -n "$transcriptome_db" ]; then
    echo "Transcriptome database: $transcriptome_db" | tee -a "$log_file"
else
    transcriptome_db="None"
fi

if [ "$orf_type" = "canonical" ] || [ "$orf_type" = "all" ]; then
    echo "ORF type: $orf_type" | tee -a "$log_file"
else
    echo "The ORF type '$orf_type' is not supported. Please specify either 'canonical' or 'all'." | tee -a "$log_file"
    exit 1
fi

if [ "$upstream_orf" = "TRUE" ] || [ "$upstream_orf" = "FALSE" ]; then
    echo "Upstream ORF detection: $upstream_orf" | tee -a "$log_file"
else
    echo "Incorrect boolean value for upstream ORF detection. Please specify either 'TRUE' or 'FALSE'." | tee -a "$log_file"
    exit 1
fi

if [ "$downstream_orf" = "TRUE" ] || [ "$downstream_orf" = "FALSE" ]; then
    echo "Downstream ORF detection: $downstream_orf" | tee -a "$log_file"
else
    echo "Incorrect boolean value for downstream ORF detection. Please provide either 'TRUE' or 'FALSE'." | tee -a "$log_file"
    exit 1
fi

if ! [ -n "$orf_length" ]; then
    orf_length=30
fi

echo "Minimum ORF length (in amino acids): $orf_length" | tee -a "$log_file"

if ! [ -n "$threads" ]; then
    threads=1
fi

echo "Number of threads: $threads" | tee -a "$log_file"

if ! [ -n "$min_tx_count" ]; then
    min_tx_count=5
fi

echo "Transcript expression cutoff: $min_tx_count" | tee -a "$log_file"

if [ -n "$tx_count_file" ] ; then
    echo "Transcript count file: $tx_count_file" | tee -a "$log_file"
else
    tx_count_file="None"
fi

fastq_id_pe=()
fastqfile_se=()
bamfiles=()
index=0

# Check if the sample directory exists
if [ -d "$sample_directory" ] ; then
    echo "Searching for $data_type files from the sample directory '$sample_directory'..." | tee -a "$log_file"

    if [ "$data_type" = "FASTQ" ]; then
        for file in "$sample_directory"/*.fastq "$sample_directory"/*.fastq.gz; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                if [[ $filename == *_R1.fastq* ]]; then
                    sample_id="${filename/_R1.fastq*/}"
                    fastq_id_pe[$index]=$sample_id
                    index=$((index+1))
                elif [[ $filename != *_R1.fastq* && $filename != *_R2.fastq* ]]; then
                    fastqfile_se[$index]=$filename
                    index=$((index+1))
                fi
            fi
        done
    elif [ "$data_type" = "BAM" ]; then
        for file in "$sample_directory"/*.bam; do
            if [ -f "$file" ]; then
                bamfiles[$index]=$(basename "$file")
                echo $file
                index=$((index+1))
            fi
        done
    fi
else
    echo "The sample directory does not exist." | tee -a "$log_file"
    exit 1
fi

# Set the working directory
cd ../../

genomeprot_dir=$(pwd)
echo "GenomeProt direcory: $genomeprot_dir"

echo $log_file
# Capture end time
start=$(date)
start_time=$(date +%s)

echo "Run started at: "$start | tee -a "$log_file"

if [[ "$data_type" == "GTF" ]]; then
    generate_database "$custom_gtf" "$output_directory" "$reference_gtf" "$organism" "$orf_length" "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$genome_fa" "$tx_count_file" "$min_tx_count" "$data_type"

elif [[ "$sequencing_platform" == "short-read" && -f "$transcriptome_db" ]]; then
    # run salmon and compile transcript and transcript_counts.txt file
    short_read_commands fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$genome_fa" "$transcriptome_db" "$reference_gtf" "$sample_directory" "$output_directory" "$organism" "$orf_length" "$vcf" "$orf_type" "$threads" "$upstream_orf" "$downstream_orf" "$ref_proteome"

    # generate database
    tx_count_file="${output_directory}salmon_transcript_counts.txt"
    generate_database "$reference_gtf" "$output_directory" "$reference_gtf" "$organism" "$orf_length" "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$genome_fa" "$tx_count_file" "$min_tx_count" "$data_type"

elif [[ "$sequencing_platform" == "long-read" ]]; then
    fastq_bam_input fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$sequencing_platform" "$genome_fa" "$sample_directory" "$output_directory" "$threads" "$reference_gtf" "$organism" "$orf_length" "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$min_tx_count" "$data_type"

    # Rename files
    mv "$output_directory/extended_annotations.gtf" "$output_directory/bambu_transcript_annotations.gtf"
    mv "$output_directory/counts_transcript.txt" "$output_directory/bambu_transcript_counts.txt"
    mv "$output_directory/counts_gene.txt" "$output_directory/bambu_gene_counts.txt"

else
    echo "Error: unsupported combination (sequencing type: '$sequencing_platform', data type: '$data_type') or missing transcriptome database." >&2
    exit 1
fi

# Capture end time
end=$(date)
end_time=$(date +%s)

echo "Run ended at: "$end | tee -a "$log_file"
duration=$((end_time - start_time))

# Convert duration to days, hours, minutes, and seconds
days=$((duration / 86400))
hours=$(( (duration % 86400) / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Runtime: $days days, $hours hours, $minutes minutes, $seconds seconds" | tee -a "$log_file"