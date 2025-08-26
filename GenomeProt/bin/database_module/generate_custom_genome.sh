#!/bin/bash

# usage() function to show the help message
function usage() {
    cat << EOF
Program: Generate custom genome

Usage: generate_custom_genome.sh [-h] [-g <genome_fasta>] [-r <reference GTF>] [-v <VCF file>] [-o <outdir>] 

[-l <ORF length (AA)>] [-j <job identifier>] [-d <sample directory>]

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

vcf_file_without_extn=$(echo $vcf_file_name| sed s'/.vcf//g')

genome_file=$( basename $genome_fa)
genome_file_without_extn=$( echo $genome_file | sed s'/.fa//g')


/opt/miniconda3/bin/conda run -n GenomeProt_env --no-capture-output python3 bin/database_module/vcfparser.py $vcf_file $outdir


/opt/miniconda3/bin/conda run -n GenomeProt_env bgzip -c $outdir$vcf_file_without_extn"_heterozygous.vcf" >$outdir$vcf_file_without_extn"_heterozygous.vcf.gz"
/opt/miniconda3/bin/conda run -n GenomeProt_env bgzip -c $outdir$vcf_file_without_extn"_homozygous.vcf" >$outdir$vcf_file_without_extn"_homozygous.vcf.gz"
/opt/miniconda3/bin/conda run -n GenomeProt_env tabix -p vcf $outdir$vcf_file_without_extn"_homozygous.vcf.gz"
/opt/miniconda3/bin/conda run -n GenomeProt_env tabix -p vcf $outdir$vcf_file_without_extn"_heterozygous.vcf.gz"
/opt/miniconda3/bin/conda run -n GenomeProt_env bcftools consensus -f $genome_fa -o $outdir$genome_file_without_extn"_hm.fa" $outdir$vcf_file_without_extn"_homozygous.vcf.gz"
/opt/miniconda3/bin/conda run -n GenomeProt_env bcftools consensus -f $outdir$genome_file_without_extn"_hm.fa" -o $outdir$genome_file_without_extn"_hm_ht.fa" $outdir$vcf_file_without_extn"_heterozygous.vcf.gz"

conda deactivate
