#!/bin/bash

# usage() function to show the help message


function genomeprot_commands() {
  fq_se=("${!1}")    # Use indirect referencing to dereference the array
  bm_fn=("${!2}")
  fastq_id_pe=("${!3}")
  platform=$4
  genome_fa=$5
  sample_dir=$6
  output_dir=$7
  threads=$8
  gtf=$9
  #echo "${fq_se[@]}"
  #echo "${bm_fn[@]}"
  #echo "${fastq_id_pe[@]}"
  echo "Sequencing platform is:\t" $platform
  index_fn=$(basename "$genome_fa")".mmi"
  
  # Apply a for loop
  if [[ $platform == "long-read-ont" ]];then
    echo minimap2 -ax splice:hq -d  $output_dir$index_fn  $genome_fa
    if [[ ${#fq_se[@]} -gt 0 ]]; then
      for file in "${fq_se[@]}"; do
        echo "Processing: $file"
        file_prefix="${file%%.fastq*}"
        echo $file_prefix
        echo eval minimap2 -t $threads -ax splice:hq --sam-hit-only --secondary=no $output_dir$index_fn $sample_dir$file  \| samtools view -bh \| samtools sort -@ $threads -o  $output_dir$file_prefix".bam"

        
      done
      
      #run bambu
      echo Rscript ./GenomeProt/bin/database_module/run_bambu.R -b $output_dir -g $gtf -f $genome_fa -o $output_dir
      
    fi
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

Usage: genomeprot_db_generation.sh [-h] [-s <sequencing platform>] [-o <organism>] [-a <annotations>] [-c <start codon>] [-l <ORF length (AA)>] [-t <data type>] [-g <reference genome fasta][-d <sample directory>] [-O <output directory>]

OPTIONS:
    -h  help options
    -s  Sequencing platform: long-read-ont/long-read-pacbio/short-read
    -o  Organism: human/mouse
    -c  Start codon: ATG, ATG+CTG
    -l  ORF length (in amino acids)
    -t  Data type: FASTQ, BAM, or GTF
    -g  Genome fasta file
    -@  Threads 
    -a  Annotaton: GENCODE GTF file 
    -d  Sample directory
    -O  Output directory
    
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
OPTS=$(getopt -o s:o:c:l:t:g:@:a:d:O:h -n "genomeprot_db_generation.sh" -- "$@")

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
        -h)
            usage
            exit 0
            ;;
        -s)
            sequencing_platform=$2
            shift 2
            ;;
        -o)
            organism=$2
            shift 2
            ;;
        -c)
            start_codon=$2
            
            shift 2
            ;;
        -l)
            orf_length=$2
            shift 2
            ;;
        -t)
            data_type=$2
            shift 2
            ;;
        -g)
            genome_fa=$2
            shift 2
            ;;
        -@)
            threads=$2
            shift 2
            ;;
        -a)
            annotations=$2
            shift 2
            ;;
        -d)
            sample_directory=$2
            shift 2
            ;;
        
        -O)
            output_directory=$2
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

# Checking input arguments are correct

#checking where actually files and folders provided by the user are correct.

if [ "$sequencing_platform" = "long-read-ont" ] || [ "$sequencing_platform" = "long-read-pacbio" ] || [ "$sequencing_platform" = "short-read" ]; then
    echo "Sequencing platform: $sequencing_platform"
else 
    echo "Incorrect value for sequencing platform. Please provide nanopore or pacbio"
    exit 1
fi

if [ "$organism" = "mouse" ] || [ "$organism" = "human" ]; then
    echo "Organism: $organism"
else 
    echo "Incorrect value for organism. Please provide mouse or human"
    exit 1
fi


if [ "$start_codon" = "ATG" ] || [ "$start_codon" = "ATG+CTG" ]; then
    echo "Start codon: $start_codon"
else 
    echo "Incorrect value for start codon. Please provide ATG or ATG+CTG"
    exit 1
fi

#check if sample directory exists

fastq_id_pe=()
fastqfile_se=()
bamfiles=()
index=0



if [ -d "$sample_directory" ] ; then
    echo "sample directory: $sample_directory exists"
    if [ "$data_type" = "FASTQ" ]; then
    
        for file in "$sample_directory"/*.fastq "$sample_directory"/*.fastq.gz; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                if [[ $filename == *_R1.fastq* ]]; then
                    sample_id="${filename/_R1.fastq*/}"
                    #echo $sample_id
                    fastq_id_pe[$index]=$sample_id
              
                    index=$((index+1))
                elif [[ $filename != *_R1.fastq* && $filename != *_R2.fastq* ]]; then
                  fastqfile_se[$index]=$filename
                  #echo $filename
                  index=$((index+1))
                fi
            fi
        done
    elif [ "$data_type" = "BAM" ]; then
    
        for file in "$sample_directory"/*.bam; do
            if [ -f "$file" ]; then
                echo $(basename "$file")
                bamfiles[$index]=$(basename "$file")
                index=$((index+1))
            fi
        done
    
    fi

else
    echo "Please provide correct sample directory name"
    exit 1
fi

#echo "${fastq_id_pe[@]}"
#echo "${bamfiles[@]}"
#echo "${fastq_id_pe[@]}"
#Directories
cd ../../..
genomeprot_dir=$(pwd)

echo "GenomeProt direcory: $genomeprot_dir"

yaml_file=$genomeprot_dir"/conda_env.yaml"

echo $yaml_file


# Source the Conda script to enable the `conda` command
conda_path="/opt/miniconda3/etc/profile.d/conda.sh"


# Name of the conda environment
genomeprot_env="GenomeProt_env"

# Check if the Conda environment exists
if conda env list | grep -q "^$genomeprot_env\s"; then
    echo "Conda environment '$genomeprot_env'exists."
    echo "Activating conda environment"
    source $conda_path
    conda activate $genomeprot_env
    #execute commands
    genomeprot_commands fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$sequencing_platform" "$genome_fa" "$sample_directory" "$output_directory" "$threads" "$annotations"
else
    echo "Conda environment '$genomeprot_env' does not exist. Creating it from $yaml_file..."
    if [ -f "$yaml_file" ]; then
        conda env create -f "$yaml_file"
        if [ $? -eq 0 ]; then
            echo "Conda environment '$genomeprot_env' created successfully."
            
            #execute commands
            
            ################
            genomeprot_commands fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$sequencing_platform" "$genome_fa" "$sample_directory" "$output_directory" "$threads" "$annotations"
            
            
            ################
        else
            echo "Failed to create Conda environment '$genomeprot_env'. Please check the YAML file." >&2
            exit 1
        fi
    else
        echo "YAML file '$yaml_file' not found. Cannot create Conda environment." >&2
        exit 1
    fi
fi






