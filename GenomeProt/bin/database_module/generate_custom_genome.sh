#!/bin/bash

function usage() {
    /usr/bin/cat << EOF
Program: Custom genome generation tool

Usage: generate_custom_genome.sh [-h] [-g <reference genome FASTA>] [-v <VCF file>] [-o <output directory>]

Generate a custom genome where genomic variants from a VCF file are injected into a reference genome.

OPTIONS
    -h              Print this help message
    -g <file>       Reference genome FASTA
    -v <file>       VCF file with variants from all samples
    -o <dir>        Output directory

=== Expected file formats ===

The reference genome FASTA (-g) must have a file extension of '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa' or '.frn'.
The VCF file (-v / --vcf) must have a file extension of '.vcf'.

EOF
}

# Get the path to the conda binary
conda_path=$(/usr/bin/which conda)
if [ -z "$conda_path" ]; then
    /usr/bin/echo 'Error: Failed to find the conda binary. Check if conda is installed or if it is within a directory specified in the $PATH variable.' >&2
    exit 1
fi

# Name of the Conda environment
genomeprot_env="GenomeProt_env"

# Check if the Conda environment exists
if ! ("$conda_path" env list | grep -q "^$genomeprot_env\s"); then
    /usr/bin/echo "Error: The Conda environment '$genomeprot_env' does not exist." >&2
    exit 1
fi

# Ensure that the script exits if an error occurs
set -e

# Check if no arguments were provided
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Parse options using getopt
OPTS=$(/usr/bin/getopt -o g:v:o:h -n "generate_custom_genome.sh" -- "$@")

if [ $? -ne 0 ]; then
    /usr/bin/echo 'Error: Failed to parse options.' >&2
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
        -g)
            genome_fa=$2
            shift 2
            ;;
        -v)
            vcf_file=$2
            shift 2
            ;;
        -o)
            outdir=$2
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            /usr/bin/echo "Invalid option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [ -z "$outdir" ]; then
    echo 'Error: An output directory must be provided.' >&2
    exit 1
fi

if [ -z "$genome_fa" ]; then
    echo 'Error: A reference genome FASTA must be provided.' >&2
    exit 1
elif ! [ -f "$genome_fa" ]; then
    echo "Error: The reference genome FASTA '$genome_fa' does not exist but is required." >&2
    exit 1
elif ! [ -s "$genome_fa" ]; then
    echo "Error: The reference genome FASTA '$genome_fa' is empty." >&2
    exit 1
elif ! ([[ "$genome_fa" == *.fasta ]] || [[ "$genome_fa" == *.fas ]] || [[ "$genome_fa" == *.fa ]] || [[ "$genome_fa" == *.fna ]] || [[ "$genome_fa" == *.ffn ]] || [[ "$genome_fa" == *.faa ]] || [[ "$genome_fa" == *.mpfa ]] || [[ "$genome_fa" == *.frn ]]); then
    echo "Error: The reference genome FASTA file '$genome_fa' does not have one of the following accepted file extensions: .fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn" >&2
    exit 1
fi

if [ -z "$vcf_file" ]; then
    echo 'Error: A VCF file must be provided.' >&2
    exit 1
elif ! [ -f "$vcf_file" ]; then
    echo "Error: The VCF file '$vcf_file' does not exist but is required." >&2
    exit 1
elif ! [ -s "$vcf_file" ]; then
    echo "Error: The VCF file '$vcf_file' is empty." >&2
    exit 1
elif ! [[ "$vcf_file" == *.vcf ]]; then
    echo "Error: The VCF file '$vcf_file' does not have the .vcf file extension." >&2
    exit 1
fi

vcf_file_name=$(/usr/bin/basename $vcf_file)
vcf_file_without_extn="${vcf_file_name%.*}"

genome_file=$(/usr/bin/basename $genome_fa)
genome_file_without_extn="${genome_file%.*}"

# Segregate homozygous and heterozygous variants in the VCF file into two separate files
"$conda_path" run -n "$genomeprot_env" python3 bin/database_module/vcfparser.py "$vcf_file" "$outdir"

# Compress the VCF files
"$conda_path" run -n "$genomeprot_env" bgzip "${outdir}/${vcf_file_without_extn}_heterozygous.vcf" -o "${outdir}/${vcf_file_without_extn}_heterozygous.vcf.gz"
"$conda_path" run -n "$genomeprot_env" bgzip "${outdir}/${vcf_file_without_extn}_homozygous.vcf" -o "${outdir}/${vcf_file_without_extn}_homozygous.vcf.gz"

# Index the VCF files
"$conda_path" run -n "$genomeprot_env" tabix -p vcf "${outdir}/${vcf_file_without_extn}_homozygous.vcf.gz"
"$conda_path" run -n "$genomeprot_env" tabix -p vcf "${outdir}/${vcf_file_without_extn}_heterozygous.vcf.gz"

# Create the consensus sequence by applying the homozygous VCF variants to the reference genome
"$conda_path" run -n "$genomeprot_env" bcftools consensus -f "$genome_fa" -o "${outdir}/${genome_file_without_extn}_hm.fa" "${outdir}/${vcf_file_without_extn}_homozygous.vcf.gz" -H 2

# Create the consensus sequence by applying the heterozygous VCF variants to the genome
if [ ! -s "${outdir}/${genome_file_without_extn}_hm_ht.fa" ]; then
    echo "File is empty"
    "$conda_path" run -n "$genomeprot_env" bcftools consensus -f "$genome_fa" -o "${outdir}/${genome_file_without_extn}_hm_ht.fa" "${outdir}/${vcf_file_without_extn}_heterozygous.vcf.gz" -H 2
else
    echo "File is not empty"
    "$conda_path" run -n "$genomeprot_env" bcftools consensus -f "${outdir}/${genome_file_without_extn}_hm.fa" -o "${outdir}/${genome_file_without_extn}_hm_ht.fa" "${outdir}/${vcf_file_without_extn}_heterozygous.vcf.gz" -H 2
fi