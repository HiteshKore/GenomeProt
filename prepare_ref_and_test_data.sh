#!/bin/sh

# Unzip the combined UniProt and OpenProt reference files

cd GenomeProt/data

for filepath in openprot_uniprotDb_c_elegans.txt openprot_uniprotDb_drosophila.txt openprot_uniprotDb_hs.txt openprot_uniprotDb_mm.txt openprot_uniprotDb_rat.txt openprot_uniprotDb_zebrafish.txt; do
    if ! [ -f "$filepath" ] && [ -f "$filepath.zip" ]; then
        unzip "$filepath.zip"
        rm "$filepath.zip"
    fi
done

cd ..

# Download the test data for the database generation module to preload

mkdir -p testdata/long_read_bam
cd testdata

for filepath in gencode_v47_sorted.gtf BRAF_mutation.vcf GRCh38_chr1_6_7.fa.gz; do
    if ! [ -f "$filepath" ]; then
        curl -O "https://genomeprot.researchsoftware.unimelb.edu.au/testdata/$filepath"
    fi
done

if ! [ -f long_read_bam/Melanoma_data_subset.bam ]; then
    curl -o long_read_bam/Melanoma_data_subset.bam https://genomeprot.researchsoftware.unimelb.edu.au/testdata/long_read_bam/Melanoma_data_subset.bam
fi

cd ../..