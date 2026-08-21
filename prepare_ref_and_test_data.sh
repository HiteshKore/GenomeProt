#!/bin/sh

# Check if the script is being run in the directory it is in

if ! [ -d GenomeProt/data ] || ! [ -d GenomeProt/testdata ]; then
    echo 'Error: Either the GenomeProt/data directory or the GenomeProt/testdata directory is not found. Please run this script in the directory it is in.' >&2
    exit 1
fi

# Unzip the combined UniProt and OpenProt reference files

cd GenomeProt/data

for filepath in openprot_uniprotDb_human.txt openprot_uniprotDb_c_elegans.txt openprot_uniprotDb_drosophila.txt openprot_uniprotDb_mouse.txt openprot_uniprotDb_rat.txt openprot_uniprotDb_zebrafish.txt; do
    if ! [ -f "$filepath" ] && [ -f "$filepath.zip" ]; then
        echo "Unzipping GenomeProt/data/$filepath.zip..."
        unzip "$filepath.zip"
        rm "$filepath.zip"
        echo 'Done!'
    fi
done

cd ..

# Unzip the test data files

cd testdata

for filepath in BRAF_mutation.vcf gencode_v47_sorted.gtf GRCh38_chr1_6_7_masked.fa.gz peptide_data.tsv; do
    if ! [ -f "$filepath" ] && [ -f "$filepath.zip" ]; then
        echo "Unzipping GenomeProt/testdata/$filepath.zip..."
        unzip "$filepath.zip"
        rm "$filepath.zip"
        echo 'Done!'
    fi
done

if ! [ -f GRCh38_chr1_6_7_masked.fa ] && [ -f GRCh38_chr1_6_7_masked.fa.gz ]; then
    echo 'Decompressing GenomeProt/testdata/GRCh38_chr1_6_7_masked.fa.gz...'
    gzip -d GRCh38_chr1_6_7_masked.fa.gz
    echo 'Done!'
fi

if ! [ -f long_read_bam/Melanoma_data_subset.bam ] && [ -f long_read_bam.zip ]; then
    echo 'Unzipping GenomeProt/testdata/long_read_bam.zip...'
    unzip long_read_bam.zip
    rm long_read_bam.zip
    echo 'Done!'
fi

cd ../..

echo 'The reference and test datasets are now unzipped and ready to use in GenomeProt.'