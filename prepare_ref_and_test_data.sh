#!/bin/sh

# Unzip the combined UniProt and OpenProt reference files

cd GenomeProt/data

for filepath in openprot_uniprotDb_human.txt openprot_uniprotDb_c_elegans.txt openprot_uniprotDb_drosophila.txt openprot_uniprotDb_mouse.txt openprot_uniprotDb_rat.txt openprot_uniprotDb_zebrafish.txt; do
    if ! [ -f "$filepath" ] && [ -f "$filepath.zip" ]; then
        unzip "$filepath.zip"
        rm "$filepath.zip"
    fi
done

cd ..

# Unzip the test data files

cd testdata

for filepath in BRAF_mutation.vcf gencode_v47_sorted.gtf GRCh38_chr1_6_7_masked.fa.gz peptide_data.tsv; do
    if ! [ -f "$filepath" ] && [ -f "$filepath.zip" ]; then
        unzip "$filepath.zip"
        rm "$filepath.zip"
    fi
done

gunzip GRCh38_chr1_6_7_masked.fa.gz

if ! [ -f long_read_bam/Melanoma_data_subset.bam ] && [ -f long_read_bam.zip ]; then
    unzip long_read_bam.zip
    rm long_read_bam.zip
fi

cd ../..