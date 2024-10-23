#!/bin/bash

# usage() function to show the help message
function usage() {
    cat << EOF
Program: Generate custom genome

Usage: genomeprot_db_generation.sh [-h] [-g <genome_fasta>] [-o <organism>] [-a <annotations>] [-c <start codon>] [-l <ORF length (AA)>] [-j <job identifier>] [-d <sample directory>]

OPTIONS:
    -h help		help options
    -g genome fasta		Reference genome FASTA file   
    -r reference GTF	Reference GTF file in GENCODE/ENSEMBL format
    -v VCF file  Varinats in VCF file format
    -o outdir		Output directory
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
OPTS=$(getopt -o g:r:v:o:h -n "generate_custom_genome.sh" -- "$@")

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
        -g)
            genome_fa=$2
            shift 2
            ;;
        -r)
            reference_gtf=$2
            
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
            echo "Invalid option: $1" >&2
            usage
            exit 1
            ;;
    esac
done


# genome_fa=/home/hiteshk/GenomeProt/PGCode/reference/GRCh38.primary_assembly_chrX.fa
# reference_gtf=/home/hiteshk/GenomeProt/PGCode/reference/gencode.v45.chr_patch_hapl_scaff.annotation_hs_chrX.gtf
# vcf_file=/home/hiteshk/GenomeProt/PGData/Mutant_peptide_Db/Cosmic_chrX.vcf.gz
# outdir=$(dirname "$vcf_file")/

vcf_file_name=$( basename $vcf_file)

# Check if the file ends with .gz
if [[ "${vcf_file}" == *.gz ]]; then
    flt_vcf_file=$( basename $vcf_file | sed s'/.vcf.gz/.vcf_norm.gz/g')
else
    echo "The file is not compressed. Compressing now..."
    bgzip -c "$vcf_file" >$outdir$vcf_file_name".gz"
    vcf_file=$outdir$vcf_file_name".gz"
    tabix -p vcf $vcf_file

    flt_vcf_file=$( echo $vcf_file_name".gz"| sed s'/.vcf.gz/.vcf_norm.gz/g')
    
    echo $vcf_file
    
    
fi

#altered genome output file name
alt_genome_file=$( basename $genome_fa | sed s'/.fa/_alt.fa/g')

#remove consider first altered allele in case of multiple altered alleles

 norm_vcf="bcftools view -v snps $vcf_file | awk 'BEGIN {FS=OFS=\"\\t\"} {split(\$5, alleles, \",\"); \$5=alleles[1]; print}' | bcftools norm -f $genome_fa -o $outdir$flt_vcf_file -c x -d all -O z"
 echo $norm_vcf

 bash -c "$norm_vcf"

#index filtered VCF file 
tabix -p vcf $outdir$flt_vcf_file
 
bcftools consensus -f $genome_fa -o $outdir$alt_genome_file $outdir$flt_vcf_file




