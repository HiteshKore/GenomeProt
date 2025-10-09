#!/bin/bash

function short_read_commands(){
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

  # generate salmon index
  # Check if both reference genome and transcriptome are gzipped
  echo "Generating salmon index"
  if [[ "$transcriptome_db" == *.gz ]] && [[ "$genome_fa" == *.gz ]]; then
      echo "Processing gzipped files"
      # Generate decoy file
      grep '^>' <(gunzip -c $genome_fa) | cut -d ' ' -f 1 > $output_dir"decoys.txt"
            
      # Remove '>' from decoy file
      sed -i -e 's/>//g' $output_dir"decoys.txt"
      
      # Concatenate reference and transcriptome to create gentrome
      cat "$transcriptome_db" "$genome_fa" > $output_dir"gentrome.fa.gz"
  
      # Run Salmon index
      conda run -n GenomeProt_env salmon index -t $output_dir"gentrome.fa.gz" -d $output_dir"decoys.txt"  -p "$threads" -i $output_dir"salmon_index" --gencode
  
  else
      echo "Processing uncompressed files"
      # Generate decoy file
      grep '^>' "$genome_fa" | cut -d " " -f 1 > $output_dir"decoys.txt"
  
      # Remove '>' from decoy file
      sed -i -e 's/>//g' $output_dir"decoys.txt"
  
      # Concatenate reference and transcriptome to create gentrome
      cat "$transcriptome_db" "$genome_fa" > $output_dir"gentrome.fa"
      # Run Salmon index
      conda run -n GenomeProt_env salmon index -t $output_dir"gentrome.fa" -d $output_dir"decoys.txt" -p "$threads" -i $output_dir"salmon_index" --gencode

fi
  

  #paired-end data
  #loop over the samples
  if [[ ${#fastq_id_pe[@]} -gt 0 ]]; then
    pe_ids=() #sample ids
    for id in "${fastq_id_pe[@]}"; do
     echo "Processing: "$id | tee -a "$log_file"
       
      r1_file_path=$(find "$sample_dir" -type f \( -name "${id}_R1.fastq" -o -name "${id}_R1.fastq.gz" \))
      r2_file_path=$(find "$sample_dir" -type f \( -name "${id}_R2.fastq" -o -name "${id}_R2.fastq.gz" \))
      pe_ids+=("$id")
      conda run -n GenomeProt_env salmon quant -i $output_dir"salmon_index" -p "$threads" -l A -1 $r1_file_path -2 $r2_file_path --validateMappings -o $output_dir$id
  
      echo $file_path
  
    done
    
    sample_string_pe="${pe_ids[*]}"
     
    echo "#generate transcript_counts file"
    conda run -n GenomeProt_env Rscript ./bin/database_module/matrix_compilation_salmon.R -q  "$output_dir" -s "$sample_string_pe" -g "$ref_gtf"
  
  fi
  
  #single-end data
  #loop over the samples
  if [[ ${#fq_se[@]} -gt 0 ]]; then
  
    
    ids=() #sample ids
    for file in "${fq_se[@]}"; do
       id=$(echo "$file" | sed -e 's/\(.*\)\.fastq\.gz/\1/' -e 's/\(.*\)\.fastq/\1/')
       echo "Processing: "$id | tee -a "$log_file"
       ids+=("$id")
       conda run -n GenomeProt_env salmon quant -i $output_dir"salmon_index" -p "$threads" -l A -r $sample_dir$file --validateMappings -o $output_dir$id
       
    done
  
  sample_string="${ids[*]}"
  #generate transcript_counts file
  conda run -n GenomeProt_env Rscript ./bin/database_module/matrix_compilation_salmon.R -q  "$output_dir" -s "$sample_string" -g "$ref_gtf"
  
  fi

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
  
  
    if [[ "$vcf_fn" == "None" ]]; then 
        #generate proteome
    
        echo "Generating proteome database" | tee -a "$log_file"
        
        echo $tx_count_file
        echo $min_tx_count
        if [[ "$tx_count_file" == "None" ]]; then
        
        conda run -n GenomeProt_env Rscript ./bin/database_module/generate_proteome.R -g $custom_gtf -r $ref_gtf -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -s $output_dir
        
        else 
        
        conda run -n GenomeProt_env Rscript ./bin/database_module/generate_proteome.R -g $custom_gtf -r $ref_gtf -c $tx_count_file -m $min_tx_count -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -s $output_dir
        
        fi
        
        if [[ "$input_type" == "GTF" ]];then
            echo "Annotating proteome database" | tee -a "$log_file"
            conda run -n GenomeProt_env --no-capture-output python ./bin/database_module/annotate_proteome.py $ref_gtf $ref_proteome $output_dir"ORFome_aa.txt" $custom_gtf $output_dir $orf_input_type $orf_length None $organism
        else
            echo "Annotating proteome database" | tee -a "$log_file"
            conda run -n GenomeProt_env --no-capture-output python ./bin/database_module/annotate_proteome.py $ref_gtf $ref_proteome $output_dir"ORFome_aa.txt" $output_dir"proteome_database_transcripts.gtf" $output_dir $orf_input_type $orf_length None $organism
        
        fi

      
    else #vcf file provided
    
        echo $tx_count_file
        echo $min_tx_count
        if [[ "$tx_count_file" == "None" ]]; then
          echo "Generating proteome database" | tee -a "$log_file"
          conda run -n GenomeProt_env Rscript ./bin/database_module/generate_proteome.R -G $ref_genome -g $custom_gtf -r $ref_gtf -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -v $vcf_fn -s $output_dir
        else
          conda run -n GenomeProt_env Rscript ./bin/database_module/generate_proteome.R -G $ref_genome -g $custom_gtf -r $ref_gtf -c $tx_count_file -m $min_tx_count -o $organism -l $orf_length -u $five_utr_orf -d $three_utr_orf -v $vcf_fn -s $output_dir
        fi
      

      echo "Annotating proteome database" | tee -a "$log_file"
      conda run -n GenomeProt_env --no-capture-output python ./bin/database_module/annotate_proteome.py $gtf $ref_proteome $output_dir"ORFome_aa.txt" $output_dir"proteome_database_transcripts.gtf" $output_dir $orf_input_type $orf_length $vcf_fn $organism

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

  echo "Sequencing platform is: " $platform | tee -a "$log_file"
  index_fn=$(basename "$genome_fa")".mmi"
  
  # Apply a for loop
  if [[ $platform == "long-read" ]];then
   
    
      if [[ ${#fq_se[@]} -gt 0 ]]; then
        #generate index
        echo "Generating minimap2 index" | tee -a "$log_file"
         # minimap2 -ax splice:hq -d  $output_dir$index_fn  $genome_fa
        conda run -n GenomeProt_env bash -c "minimap2 -ax splice:hq -d  $output_dir$index_fn  $genome_fa"
        echo "Indexing complete" | tee -a "$log_file"
        echo "Genome alignment started" | tee -a "$log_file"
        for file in "${fq_se[@]}"; do
          file_prefix="${file%%.fastq*}"
          echo "Processing:\t"$file_prefix | tee -a "$log_file"
          conda run -n GenomeProt_env bash -c "minimap2 -t $threads -ax splice:hq --sam-hit-only --secondary=no $output_dir$index_fn $sample_dir$file | samtools view -bh | samtools sort -@ $threads -o  $output_dir$file_prefix.bam"
          echo "Transcript identification/quantification using bambu" | tee -a "$log_file"
          
        done
        echo "Genome alignment completed" | tee -a "$log_file"
        conda run -n GenomeProt_env Rscript ./bin/database_module/run_bambu.R -b $output_dir -g $ref_gtf -s $organism -o $output_dir
      fi #fastq condition close
     
    
      #run bambu
      if [[ ${#bm_fn[@]} -gt 0 ]]; then
      
        echo "Transcript identification/quantification using bambu" | tee -a "$log_file"
        conda run -n GenomeProt_env Rscript ./bin/database_module/run_bambu.R -b $sample_dir -g $ref_gtf -s $organism -o $output_dir
      
      fi
      
      
      transcript_count_fn=$output_directory"counts_transcript.txt"
      custom_gtf=$output_directory"extended_annotations.gtf"
      
      echo "****"$transcript_count_fn
      echo "****"$transcript_expr_cutoff
   
      #protome database generation
      generate_database "$custom_gtf" "$output_dir" "$ref_gtf" "$organism" "$orf_length" "$vcf_fn" "$orf_input_type" "$five_utr_orf" "$three_utr_orf" "$ref_proteome" "$genome_fa" "$transcript_count_fn" "$transcript_expr_cutoff" "$input_type"

  
  fi #platform loop closed
  
}




function usage() {
    cat << EOF
Program: GenomeProt Database Module Commandline
Version: 0.0.1v

 ▗▄▄▖▗▄▄▄▖▗▖  ▗▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▖ ▗▄▄▖  ▗▄▖▗▄▄▄▖
▐▌   ▐▌   ▐▛▚▖▐▌▐▌ ▐▌▐▛▚▞▜▌▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌ █  
▐▌▝▜▌▐▛▀▀▘▐▌ ▝▜▌▐▌ ▐▌▐▌  ▐▌▐▛▀▀▘▐▛▀▘ ▐▛▀▚▖▐▌ ▐▌ █  
▝▚▄▞▘▐▙▄▄▖▐▌  ▐▌▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▌   ▐▌ ▐▌▝▚▄▞▘ █  

Usage: genomeprot_db_generation.sh [-h] [-s <sequencing platform>] [-o <organism>] [-a <reference GTF>] [-c <start codon>] [-l <ORF length (AA)>] [-t <data type>] [-g <reference genome fasta>] [-p <threads>] [-a <reference GTF>] [-G <custom GTF>] [-v <VCF file>] [-T <ORF type>] [-r <transcriptome database>] [-U <upstream ORFs>] [-D <downstream ORFs>] [-m <minimum transcript count>] [-C <transcript count file>] [-d <sample directory>] [-O <output directory>]

OPTIONS:
    -h, --help                  Display this help and exit
    -s, --sequencing-type       Sequencing platform: long-read/short-read
    -o, --organism              Organism: HUMAN, MOUSE, RAT, CAEEL, DROME, DANRE
    -c, --start-codon           Start codon: ATG, ATG+CTG
    -l, --orf-length            ORF length in amino acids (numeric)
    -t, --data-type             Data type: FASTQ, BAM, or GTF
    -g, --genome                Reference genome FASTA file (required for FASTQ/BAM input)
    -p, --threads               Number of threads (numeric)
    -a, --reference-gtf         Reference GTF file (only supports GENCODE annotations)
    -G, --custom-gtf            Custom GTF file
    -v, --vcf                   VCF file with variants from all samples
    -T, --orf-type              ORF type: canonical/all
    -r, --transcriptome-db      Transcriptome database (required for short-read RNA-Seq input)
    -U, --upstream-orf          Upstream (5'UTR) ORFs  : TRUE/FALSE (boolean)
    -D, --downstream-orf        Downstream (3'UTR) ORFs: TRUE/FALSE (boolean)
    -m, --min-tx-count          Minimum transcript expression cutoff (numeric,default: 5)
    -C, --tx-count-file         Transcript count file (optional)
    -d, --sample-dir            Sample directory (required for FASTQ and BAM input)
    -O, --output-dir            Output directory

EOF
    exit 0
    
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

OPTS=$(getopt -o s:o:c:l:t:g:p:a:G:v:T:r:U:D:m:C:d:O:h \
              -l sequencing-type:,organism:,start-codon:,orf-length:,data-type:,genome:,threads:,reference-gtf:,custom-gtf:,vcf:,orf-type:,transcriptome-db:,upstream-orf:,downstream-orf:,min-tx-count:,tx-count-file:,sample-dir:,output-dir:,help \
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
        -c | --start-codon )        start_codon="$2"; shift 2 ;;
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


# Checking input arguments are correct

if [[ ! "$output_directory" =~ /$ ]]; then
    output_directory="$output_directory/"
fi

if [[ ! "$sample_directory" =~ /$ ]]; then
    sample_directory="$sample_directory/"
fi

#log file
log_file=$output_directory"logfile.txt"

echo $log_file

if [ "$sequencing_platform" = "long-read" ] || [ "$sequencing_platform" = "short-read" ]; then
    echo "Sequencing platform: $sequencing_platform" | tee -a "$log_file"
else 
    echo "Incorrect value for sequencing platform. Please provide long-read or short-read" | tee -a "$log_file"
    exit 1
fi


if [ "$organism" = "MOUSE" ] || [ "$organism" = "HUMAN" ] || [ "$organism" = "CAEEL" ] || [ "$organism" = "DROME" ] || [ "$organism" = "RAT" ] || [ "$organism" = "DANRE" ]; then
    
    #proteome database
  echo "Organism: $organism" | tee -a "$log_file"
  
    # set reference protein database per organism 
    if [[ "$organism" == "HUMAN" ]]; then
      ref_proteome="./data/openprot_uniprotDb_hs.txt"
    elif [[ "$organism" == "MOUSE" ]]; then
        ref_proteome="./data/openprot_uniprotDb_mm.txt"
    elif [[ "$organism" == "CAEEL" ]]; then
        ref_proteome="./data/openprot_uniprotDb_c_elegans.txt"
    elif [[ "$organism" == "DROME" ]]; then
        ref_proteome="./data/openprot_uniprotDb_drosophila.txt"
    elif [[ "$organism" == "rat" ]]; then
        ref_proteome="./data/openprot_uniprotDb_rat.txt"
    elif [[ "$organism" == "DANRE" ]]; then
        ref_proteome="./data/openprot_uniprotDb_zebrafish.txt"
    fi #organism check
  
else 
    echo "Incorrect value for organism. Please provide: HUMAN, MOUSE, CAEEL, DROME, DANRE" | tee -a "$log_file"
    exit 1
fi


if [ "$start_codon" = "ATG" ] || [ "$start_codon" = "ATG+CTG" ]; then
    echo "Start codon: $start_codon"
else 
    echo "Incorrect value for start codon. Please provide ATG or ATG+CTG" | tee -a "$log_file"
    exit 1
fi

if [ -n "$vcf" ]; then
    echo "VCF file provided: $vcf" | tee -a "$log_file"
else
    vcf="None"
fi

if [ -n "$transcriptome_db" ]; then
    echo "Transcriptome database provided: $transcriptome_db" | tee -a "$log_file"
else
    transcriptome_db="None"
fi

if [ "$orf_type" = "canonical" ] || [ "$orf_type" = "all" ]; then
    echo "ORF type: $orf_type" | tee -a "$log_file"
else 
    echo "Incorrect value for ORF type. Please provide canonical or all" | tee -a "$log_file"
    exit 1
fi

if [ "$upstream_orf" = "TRUE" ] || [ "$upstream_orf" = "FALSE" ]; then
    echo "Upstream ORF detection: $upstream_orf" | tee -a "$log_file"
else 
    echo "Incorrect value for upstream ORF detection. Please provide either TRUE or FALSE" | tee -a "$log_file"
    exit 1
fi


if [ "$downstream_orf" = "TRUE" ] || [ "$downstream_orf" = "FALSE" ]; then
    echo "Upstream ORF detection: $downstream_orf" | tee -a "$log_file"
else 
    echo "Incorrect value for downstream ORF detection. Please provide either TRUE or FALSE" | tee -a "$log_file"
    exit 1
fi


if [ -n "$min_tx_count" ]; then
    echo "Transcript expression cutoff : $tx_count_file"
else
  min_tx_count=5
fi


if [ -n "$tx_count_file" ] ; then
    echo "Transcript count file : $tx_count_file"
else 
    tx_count_file="None"
fi


if [ -n "$vcf" ]; then
    echo "VCF file provided: $vcf" | tee -a "$log_file"
else
    vcf="None"
fi

#check if sample directory exists

fastq_id_pe=()
fastqfile_se=()
bamfiles=()
index=0

if [ -d "$sample_directory" ] ; then
    echo "sample directory: $sample_directory exists" | tee -a "$log_file"
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
                index=$((index+1))
            fi
        done
    
    fi

else
    echo "Please provide correct sample directory name" | tee -a "$log_file"
    exit 1
fi

#set working directory
cd ../../../
yml_dir=$(pwd)



yaml_file=$yml_dir"/conda_env.yaml" #fix this

install_file=$yml_dir"/install.R"



echo $yaml_file 

#genomeprot directory

cd "./GenomeProt"

genomeprot_dir=$(pwd)
echo "GenomeProt direcory: $genomeprot_dir"

# Name of the conda environment
genomeprot_env="GenomeProt_env"



echo $log_file
# Capture end time
start=$(date)
start_time=$(date +%s)

echo "Run started at: "$start | tee -a "$log_file"

# Check if the Conda environment exists
if conda env list | grep -q "^$genomeprot_env\s"; then
    echo "Conda environment '$genomeprot_env'exists." 
   
    #short read data #######
    
    if [[ "$sequencing_platform" == "short-read"  && -f "$transcriptome_db" ]]; then 
        
        #run salmon and compile transcript and transcript_counts.txt file
        short_read_commands fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$genome_fa" "$transcriptome_db" "$reference_gtf" "$sample_directory" "$output_directory" "$organism" "$orf_length" "$vcf" "$orf_type" "$threads" "$upstream_orf" "$downstream_orf" "$ref_proteome"
        
        #generate database
        tx_count_file=$output_directory"counts_transcript.txt"
        generate_database "$reference_gtf" "$output_directory" "$reference_gtf" "$organism" "$orf_length"  "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$genome_fa" "$tx_count_file" "$min_tx_count" "$data_type"
        
      
    fi
       
    
    #execute commands
    if [[ "$data_type" = "GTF" ]]; then
      generate_database "$custom_gtf" "$output_directory" "$reference_gtf" "$organism" "$orf_length"  "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$genome_fa" "$tx_count_file" "$min_tx_count" "$data_type"
    else
      fastq_bam_input fastqfile_se[@] bamfiles[@] fastq_id_pe[@] "$sequencing_platform" "$genome_fa" "$sample_directory" "$output_directory" "$threads" "$reference_gtf" "$organism" "$orf_length" "$vcf" "$orf_type" "$upstream_orf" "$downstream_orf" "$ref_proteome" "$min_tx_count" "$data_type"
    fi
    
    
else
    echo "Conda environment '$genomeprot_env' does not exist. Creating it from $yaml_file..."
    if [ -f "$yaml_file" ]; then
        
        #create conda environment and install packages
        
        echo Rscript $install_file $yaml_file
        
        
    
    else
        echo "YAML file '$yaml_file' not found. Cannot create Conda environment." >&2
        exit 1
    fi #yml file condition
    
fi #environment check

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
