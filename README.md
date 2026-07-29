<img src="https://github.com/HiteshKore/GenomeProt/blob/main_beta/logo.png"/>

# GenomeProt: an integrated proteogenomics analysis platform for long-read RNA-Seq datasets

## Contents

- [Overview](#overview)
- [Installation](#installation)
  - [Option 1 (recommended): Access GenomeProt online](#option-1-recommended-access-genomeprot-online)
  - [Option 2: Run the shiny application with Docker](#option-2-run-the-shiny-application-with-docker-work-in-progress)
  - [Option 3: Locally install the shiny application](#option-3-locally-install-the-shiny-application)
- [General usage](#general-usage)
  - [Database generation](#1-database-generation)
  - [Analyse MS proteomics](#2-analyse-ms-proteomics)
  - [Integration](#3-integration)
  - [Visualisation](#4-visualisation)
- [Detailed input and output descriptions](#detailed-input-and-output-descriptions)
  - [Database generation](#1-database-generation-1)
  - [Analyse MS proteomics](#2-analyse-ms-proteomics-1)
  - [Integration](#3-integration-1)
  - [Visualisation](#4-visualisation-1)

## Overview

GenomeProt is a tool providing a web-based graphical user interface for integrative proteogenomics analyses. It generates a variant-aware, sample-specific proteome database from bulk RNA-Seq data and a VCF file, then uses it to search mass spectrometry data and identify peptides. Next, it integrates transcriptomics and proteomics results, where the identified peptides are mapped to spliced transcript coordinates. The proteogenomics integration outputs an HTML summary report on the transcripts, peptides and proteins identified, as well as BED12 files for visualisation in the UCSC Genome Browser and a combined GTF file for interactive visualisation within the GenomeProt application.

## Installation

GenomeProt can both be run as an interactive web application and run entirely from the command-line. This is because the web application is essentially a graphical wrapper around command-line scripts that perform integrative proteogenomics analyses.

To use GenomeProt now, follow [Option 1](#option-1-recommended-access-genomeprot-online).

To install GenomeProt locally or on an HPC (high-performance computing) system, follow [Option 2 (run with Docker)](#option-2-run-the-shiny-application-with-docker-work-in-progress) or [Option 3](#option-3-locally-install-the-shiny-application).

### Option 1 (recommended): Access GenomeProt online
Click on the following link to access GenomeProt in your current browser tab: https://genomeprot.researchsoftware.unimelb.edu.au/

Note: To ensure fair use of the resources, the public GenomeProt server does not perform read mapping from FASTQs (BAM or GTF input only) and excludes the "Analyse MS proteomics" step. Users wishing to perform these steps should utilise a local version of GenomeProt or perform these steps externally.

### Option 2: Run the shiny application with Docker (work in progress)
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```bash
docker run --rm -p 3838:3838 josieg/genomeprot:v1
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it's stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.

### Option 3: Locally install the shiny application

The application has substantial dependencies that we have provided as conda environment files.

#### Clone the GenomeProt GitHub repository

```bash
git clone https://github.com/HiteshKore/GenomeProt.git
```

#### Set up the GenomeProt conda environment

Install conda by following [this guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Two conda environments are provided: one contains the dependencies for running FragPipe and the GenomeProt proteomics module, while the other lacks them. Build the conda environment that suits your needs best:
```bash
cd GenomeProt
conda env create -f conda_env.yaml             # Run this line if you wish to use the proteomics module
conda env create -f conda_env_no_fragpipe.yaml # Run this line if you DO NOT wish to use the proteomics module
```

If you choose NOT to use the proteomics module, skip the step below and head to [the next step](#prepare-the-reference-and-test-datasets).

#### Install FragPipe and configure the proteomics module

The proteomics module requires a complete FragPipe v23.1 installation, which can be done by running an interactive Python script provided in this repository:

```bash
python3 fragpipe_installer.py                # Run this line if you wish to install FragPipe in the current working directory
python3 fragpipe_installer.py path/to/folder # Run this line if you wish to install FragPipe at a specific path
```

When running the script, type 'Y' and press enter to answer yes to prompts on installing FragPipe and additional necessary tools. Afterwards, follow the instructions on obtaining a 6-digit token for installing the tools:

<details>
<summary>How to obtain a 6-digit token for installing FragPipe tools <b><u>(click to expand)</u></b></summary>

1. Head to https://msfragger-upgrader.nesvilab.org/upgrader/.
2. Enter your first name, last name, academic email address and academic institution.
3. Check all tickboxes for the academic license, license agreement and SDK library distribution conditions.
4. Click on the 'Download' button (NOT the 'Get a license key...' button).
5. Wait for an email from no-reply@fragpipe.info at the email address you have specified. It should contain a download link with the following format:
  - `https://msfragger-upgrader.nesvilab.org/upgrader/download.php?token=<6-digit token>&download=<version>%24zip`
6. Copy the `<6-digit token>` from the download link and paste it into the terminal.
</details>

After the installation is complete, the following line in the `server.R` script needs to be modified if FragPipe was not extracted into `/home/user/Desktop/GenomeProt`:

```R
    " --fragpipe_path ", shQuote("/home/user/Desktop/GenomeProt/fragpipe-23.1/"),
```

For instance, if FragPipe was extracted into `/home/user/a/b/c/d/`, change `/home/user/Desktop/GenomeProt/fragpipe-23.1/` to `/home/user/a/b/c/d/fragpipe-23.1/`.

<details>
<summary>Example run of the FragPipe installation script <b><u>(click to expand)</u></b></summary>

    Extraction directory: /home/user/Desktop/GenomeProt
    Would you like to install FragPipe ([Y]/n)? Y
    ...
    Downloading FragPipe v23.1 from https://github.com/Nesvilab/FragPipe/releases/download/23.1/FragPipe-23.1-linux.zip...
    ...
    Done!
    Would you like to install MSFragger, IonQuant and diaTracer (required additional tools to run the FragPipe analysis) ([Y]/n)? Y
    ...
    To install the FragPipe tools, first head to https://msfragger-upgrader.nesvilab.org/upgrader/.
    Then, enter your first name, last name, academic email address and academic institution, and check all tickboxes for the academic license, license agreement and SDK library distribution conditions.
    Next, click on the 'Download' button (NOT the 'Get a license key...' button), and wait for an email from no-reply@fragpipe.info at the email address you have specified. It should contain a download link with the following format:
    https://msfragger-upgrader.nesvilab.org/upgrader/download.php?token=<6-digit token>&download=<version>%24zip
    The 6-digit token can be used to download MSFragger, IonQuant and diaTracer.
    Please enter your token (6 digits): 123456
    Using token '123456'...
    ...
    MSFragger downloaded!
    ...
    IonQuant downloaded!
    ...
    diaTracer downloaded!
    ...
    FragPipe installation complete.
    FragPipe extracted into /home/user/Desktop/GenomeProt.
    FragPipe tools extracted into /home/user/Desktop/GenomeProt/fragpipe-23.1/tools.
</details>

#### Prepare the reference and test datasets

Run the following command to unzip the combined UniProt and OpenProt reference files in the `GenomeProt/data` directory and download the test data for the database generation module to preload:

```bash
bash prepare_ref_and_test_data.sh
```

#### Ensuring sufficient space for storing temporary files (if required)

<details>
<summary>Instructions <b><u>(click to expand)</u></b></summary>

Files selected through the GenomeProt website are copied into a folder reserved for temporary files (usually `/tmp`). If the filesystem storing this folder is too small, selected files cannot be fully copied into the folder, which causes errors during data processing. This problem can be encountered when uploading large files, including mass spectrometry data files for FragPipe to process.

There are two methods to resolve this issue:

1. Store temporary files in a filesystem with sufficient space by setting the `TMPDIR` environment variable (does NOT require `sudo` privileges)

In this example, running `df -h` in a terminal produces the following output:

```
Filesystem      Size  Used Avail Use% Mounted on
udev            3.9G     0  3.9G   0% /dev
tmpfs           795M  524K  794M   1% /run
/dev/sda2       124G   50G   68G  43% /
tmpfs           3.9G     0  3.9G   0% /dev/shm
tmpfs           3.9G  8.0K  3.9G   1% /tmp
/dev/sda1       511M  4.7M  507M   1% /boot/efi
tmpfs           795M   36K  794M   1% /run/user/1000
```

Here, `/tmp` is its own filesystem and has a size of 3.9 GB, but that might be too small for storing mass spectrometry data files, whereas the `/dev/sda2` disk filesystem has a size of 124 GB and contains 68 GB of free space. So, the `/dev/sda2` filesystem should have enough space for storing temporary files.

Assuming the user can create files and folders in the `/dev/sda2` filesystem, they can run GenomeProt as follows, where `/path/to/temp/folder` resides in that filesystem and overrides the default location where temporary files are written to:

```bash
conda activate GenomeProt_env
mkdir -p /path/to/temp/folder
TMPDIR=/path/to/temp/folder Rscript -e "shiny::runApp('path/to/app/GenomeProt/', host='0.0.0.0', port=3838)"
```

After stopping GenomeProt by pressing Ctrl+C in the terminal where it is running, all temporary files written into `/path/to/temp/folder` should be automatically removed.

2. Increase the size of the filesystem for storing temporary files (requires `sudo` privileges; might be unsuitable for HPC systems)

Instead of changing where temporary files are written to, the filesystem storing these files can be resized to contain sufficient space. Following the example from above, assuming the user has `sudo` privileges, they can run the following command to resize the `/tmp` filesystem to 10 GB:

```bash
sudo mount -o remount,size=10G /tmp
```

This method might not work on HPC systems as such privileges are typically not provided to HPC users. Also, the recommended size of the filesystem depends on the sizes of files you wish to process through the GenomeProt website.
</details>

#### Running GenomeProt

Activate the conda environment and run the GenomeProt shiny app from the command line:
```bash
conda activate GenomeProt_env
Rscript -e "shiny::runApp('path/to/app/GenomeProt/', host='0.0.0.0', port=3838)"
```

The GenomeProt website should then be hosted locally on http://0.0.0.0:3838/ and accessible with a web browser.

To stop the web application, press Ctrl+C in the terminal tab where it is running.

## General usage

GenomeProt is an integrated proteogenomics platform with four modules: 1) database generation, 2) proteomics, 3) integration, and 4) visualisation.

<img src="https://github.com/HiteshKore/GenomeProt/blob/main_beta/GenomeProt/www/images/workflow.png"/>

### 1. Database generation

The first module generates a custom proteome database to perform proteomics searches. The module accepts RNA sequencing FASTQ files, BAM files or GTF annotation files from both short-read and long-read sequencing platforms (Illumina, Oxford Nanopore and PacBio). The main output from this module is a FASTA file with candidate protein sequences and a metadata file with details of each candidate protein.

For long-read data, transcript discovery is performed with Bambu. For short-read data, no discovery steps are performed, transcripts are instead directly quantified based on the reference transcriptome using Salmon. GenomeProt currently supports open reading frame (ORF) identification and database generation for six model organisms: fruit fly, roundworm, zebrafish, rat, mouse, and human. Users can specify an option to include short upstream ORF (uORF) and downstream ORF (dORF) protein sequences >10 amino acids (AA). Protein sequences are generated based on a user defined minimum length set to >30 AA by default. Users can also optionally provide a VCF file to incorporate single nucleotide variants (SNVs) into the genome to generate variant protein sequences.

#### Inputs:

One of the following:

 a) RNA sequencing FASTQ files, or
 b) BAM files, or
 c) a GTF annotation file

Supported sequencing platforms include short-read Illumina, long-read Oxford Nanopore and long-read PacBio.

Optional input:

- VCF file to incorporate single nucleotide variants (SNVs)

#### Outputs:

- FASTA file with candidate protein sequences
- Metadata file detailing each candidate protein

#### Features:

- Long-read transcript discovery with Bambu
- Short-read transcript quantification with Salmon
- Supports ORF identification for six model organisms: fruit fly, roundworm, zebrafish, rat, mouse, and human
- Includes optional uORF and dORF protein sequences

### 2. Analyse MS proteomics

The proteomics module uses FragPipe to perform a proteomics search with the custom database generated from the previous module. A simple interface is provided for the user to upload the custom database and mass spectrometry data files, and to configure basic FragPipe settings, which includes the proteases used, whether to add contaminants into the database and to perform peptide quantification, the number of CPU threads and the maximum amount of memory to use. Once FragPipe has finished running, the module outputs the peptide search results and peptide counts (if peptide quantification were done). Users wishing to configure FragPipe beyond the level offered by GenomeProt or to use a proteomics pipeline other than FragPipe can perform their proteomics searches externally.

#### Outputs:

- Discovered peptides from FragPipe
- Peptide counts from FragPipe (if peptide quantification was enabled)

### 3. Integration

This module integrates proteomics and transcriptomics data. Peptides are associated back to transcript isoforms and mapped to spliced genomic coordinates for downstream visualisation. This generates BED12 file of transcripts, ORFs and peptides for visualisation in the UCSC Genome Browser and a combined GTF file for visualisation within the app. An HTML report is also created that provides a summary of identified known and novel transcripts, uniquely mapping peptides, and known and novel ORFs.

#### Key outputs:

- Peptide-to-transcript mappings with spliced genomic coordinates.
- BED12 files for visualisation in the UCSC Genome Browser.
- Combined GTF file for visualisation within GenomeProt.
- HTML report summarising identified transcripts, peptides and ORFs.

### 4. Visualisation

This module uses [IsoVis](https://isomix.org/isovis/) to create an interactive visualisation that shows peptide mapping plots along transcript isoforms with quantitative peptide intensities and transcript expression data. This allows users to visualise transcript and peptide abundance across different experimental conditions. To visualise their data in IsoVis, users need to select the `combined_annotations.gtf` file generated in the integration module, and optionally transcript counts from Module 1 and peptide intensities from external proteomics analyses. Once the selected data files are processed by IsoVis, users can apply the gene filtering options provided by the webserver to quickly search for features of interest. IsoVis is highly configurable and can be used within GenomeProt or externally.

#### Features:

- Requires `combined_annotations.gtf` or `transcripts_and_ORFs_for_isovis.gtf` (from Module 3).
- Optionally input transcript isoform counts and peptide intensities.
- Allows export of plots as PDFs, SVGs, PNGs and JPEGs.

## Detailed input and output descriptions

### 1. Database generation

<table>
 <tr>
  <th>Input</th>
  <th>File Type</th>
  <th>Required?</th>
  <th>Description</th>
 </tr>
 <tr>
  <td rowspan="3">Long-read or short-read sequencing data</td>
  <td>FASTQ(s)</td>
  <td rowspan="3">One of these is required</td>
  <td>RNA sequencing reads</td>
 </tr>
 <tr>
  <td>BAM(s)</td>
  <td>Genome aligned reads</td>
 </tr>
 <tr>
  <td>GTF</td>
  <td>Assembled transcripts (from Bambu) (long-read only)</td>
 </tr>
 <tr>
  <td>Single nucleotide variants</td>
  <td>VCF</td>
  <td>No</td>
  <td>Single nucleotide variant information</td>
 </tr>
 <tr>
  <td>Transcript counts</td>
  <td>TXT/CSV</td>
  <td>No</td>
  <td>Transcript expression counts used to optionally filter lowly expressed transcripts from the database</td>
 </tr>
 <tr>
  <td>Reference annotations</td>
  <td>GTF</td>
  <td>Yes</td>
  <td>ENSEMBL or Gencode annotation</td>
 </tr>
 <tr>
  <td>Reference genome</td>
  <td>FASTA</td>
  <td>FASTQ input only</td>
  <td>Genome sequences</td>
 </tr>
 <tr>
  <td>Reference transcriptome</td>
  <td>FASTA</td>
  <td>Short-read FASTQ input only</td>
  <td>Transcriptome sequences</td>
 </tr>
</table>


| Output | File | File Type | Description |
|-------------------------------------|------------------------------------|------------|---------------------------------------------------------------------------------------------|
| Database | proteome_database.fasta | FASTA | Amino acid sequences of all ORFs in the data |
| Database metadata | proteome_database_metadata.txt | TXT | Information on each ORF in the data |
| Database transcripts (long-read data only) | proteome_database_transcripts.gtf | GTF | Annotations of transcripts used to generate the database |
| Transcript counts (FASTQ or BAM input only) | transcript_counts.txt | TXT | Transcript count file with samples as columns and transcripts as rows |
| Bambu transcript class codes (long-read FASTQ or BAM input only) | novel_transcript_classes.csv | CSV | Bambu transcript classification (see Bambu documentation) |
| GFFcompare transcript class codes (long-read FASTQ or BAM input only) | gffcompare.tmap.txt | TXT | GFFcompare transcript classification (see GFFcompare documentation) |
| Genome aligned reads (long-read FASTQ only) | sample1.bam, sample2.bam | BAM | Genome aligned reads. Only output if the input was FASTQ reads |


Proteome FASTA examples and header formats:
```
>protein_accession|CO=genomic_coordinates GA=gene_accession GN=gene_name TA=transcript_accession
MCGNNMSAPMPAVVPAARKATAAVIFLHGLGDTGHGWAEAFAGIKSPHIKYICPHAPVMPVTLNMNMAMPSWFDIVGLSPDSQEDESGIKQAAETVKALIDQEVKNGIPSNRIILGGFSQGPINSANRDISVLQCHGDCDPLVPLMFGSLTVERLKALINPANVTFKIYEGMMHSSCQQEMMDVKHFIDKLLPPID
>P10711|CO=chr1:4928137-4966584 GA=ENSMUSG00000033813.16 GN=Tcea1 TA=ENSMUST00000081551.14
MEDEVVRIAKKMDKMVQKKNAAGALDLLKELKNIPMTLELLQSTRIGMSVNALRKQSTDEEVTSLAKSLIKSWKKLLDGPSTDKDPEEKKKEPAISSQNSPEAREESSSSSNVSSRKDETNARDTYVSSFPRAPSTSDSVRLKCREMLAAALRTGDDYVAIGADEEELGSQIEEAIYQEIRNTDMKYKNRVRSRISNLKDAKNPNLRKNVLCGNIPPDLFARMTAEEMASDELKEMRKNLTKEAIREHQMAKTGGTQTDLFTCGKCKKKNCTYTQVQTRSADEPMTTFVVCNECGNRWKFC
>ORF_3|CO=chr2:53029193-53081430 GA=ENSMUSG00000061136.17 GN=Prpf40a TA=ENSMUST00000209364.3
MQATPSEAGGESPQSCLSVSRSDWTVGKPVSLLAPLIPPRSSGQPLPFGPGGRQPLRSLLVGMCSGSGRRRSSLSPTMRPGTGAERGGLMMGHPGMHYAPMGMHPMGQRANMPPVPHGMMPQMMPPMGG
```
**Note:** Unannotated ORFs are denoted by "ORF_" and variant proteins by "mORF_" followed by a unique number. UniProt or RefSeq accessions are retained for annotated proteins.

#### Open reading frame (ORF) category definitions:

| Type | Definition |
|-|-|
| CDS | Annotated in UniProt or RefSeq |
| 5UTR | Coordinates are within the 5' UTR region of an mRNA transcript |
| 3UTR | Coordinates are within the 3' UTR region of an mRNA transcript |
| 5UTR:CDS | Start site is within the 5' UTR region and stop site is within the CDS region of an mRNA transcript |
| gene_overlap | Encoded by a transcript that overlaps a region with annotated protein-coding genes |
| intergenic | Encoded by a transcript that does not overlap a region with annotated protein-coding genes |


### 2. Analyse MS proteomics

We recommend installing and running [FragPipe](https://github.com/Nesvilab/FragPipe) for analysing mass spectrometry-based proteomics data.

| Input | File Type | Required? | Description |
|-|-|-|-|
| Mass spec data | mzML, mzXML, mgf, mzBIN, Thermo RAW, Bruker PASEF .d | Yes | Mass spec files |
| Database (proteome_database.fasta) | FASTA | Yes | Generated in Module 1. Amino acid sequences of all ORFs in the data |

The output includes `peptide.tsv` and `report.pr_matrix.tsv` (the latter included only if peptide quantification were enabled by the user).


### 3. Integration

| Input | File Type | Required? | Description |
|-|-|-|-|
| Proteomics peptide data | TSV/TXT | Yes | Peptide results. Typically, 'peptides.txt', 'peptide.tsv' or 'report.pr_matrix.tsv' |
| Database (proteome_database.fasta) | FASTA | Yes | Generated in Module 1. Amino acid sequences of all ORFs in the data |
| Database metadata (proteome_database_metadata.txt) | TXT | Yes | Generated in Module 1. Information on each ORF in the data |
| Database transcripts (proteome_database_transcripts.gtf) | GTF | Yes | Generated in Module 1. Annotations of transcripts used to generate the database |


| Output | File | File Type | Description |
|-|-|-|-|
| Peptide information | peptide_info.csv | CSV | Main results file with peptide mapping data |
| Report | summary_report.html | HTML | Summary report |
| Combined annotation data | combined_annotations.gtf | GTF | Annotations of peptides, ORFs, and transcripts for visualisation |
| Transcript and ORF annotation data | transcripts_and_ORFs_for_isovis.gtf | GTF | Annotations of ORFs and transcripts for visualisation with IsoVis |
| Peptide coordinates | peptides.bed12 | BED12 | Peptide spliced genomic coordinates |
| ORF coordinates | ORFs.bed12 | BED12 | ORF spliced genomic coordinates |
| Transcript coordinates | transcripts.bed12 | BED12 | Transcript spliced genomic coordinates |


#### Description of `peptides_info.csv` output:

| Column Name | Description | Class |
|-|-|-|
| peptide | Peptide sequence | character |
| accession | Protein accession | character |
| PID | protein_accession\|CO=genomic_coordinates (included for compatibility with FASTA header) | character |
| transcript_id | ENSEMBL or novel transcript ID | character |
| gene_id | ENSEMBL gene ID | character |
| gene_name | Gene name/symbol | character |
| strand | Strand (+ or -) | character |
| number_exons | Number of exons spanned by the peptide | integer |
| transcript_length | Transcript length (nt) | integer |
| transcript_biotype | Transcript biotype from Gencode | character |
| simplified_biotype | Simplified transcript biotype | character |
| protein_length | Protein length (AA) | integer |
| orf_genomic_coordinates | Genomic coordinates of ORF | numeric |
| orf_type | Annotated, unannotated, or variant protein | character |
| localisation | ORF location in the genome (see ORF category definitions) | character |
| uniprot_status | Review status in UniProt (reviewed/unreviewed) | character |
| openprot_id | OpenProt ID if present | character |
| molecular_weight(kDA) | Molecular weight of protein (KDa) | numeric |
| isoelectric_point | Isoelectric point of ORF calculated using pKa scale EMBOSS (Rice et al., 2000) | numeric |
| hydrophobicity | Hydrophobicity profile of ORF calculated using Kyte-Doolittle scale (Kyte et al., 1982) | numeric |
| aliphatic_index | Aliphatic index of ORF (Ikai 1980) | numeric |
| longest_orf_in_transcript | Longest ORF in the transcript (longest within CDS regions for known proteins) | true/false |
| peptide_ids_gene | Is peptide uniquely mapped to gene? | true/false |
| peptide_ids_orf | Is peptide uniquely mapped to ORF? | true/false |
| peptide_ids_transcript | Is peptide uniquely mapped to transcript? | true/false |
| shared_novel_protein_peptide | Is peptide shared with other novel proteins? | true/false |
| orf_identified | Is ORF identified with unique peptide evidence? | true/false |
| gene_identified | Is gene identified with unique peptide evidence? | true/false |
| transcript_identified | Is transcript identified with unique peptide evidence? | true/false |

### 4. Visualisation

#### Peptides Visualisation using IsoVis

| Input | File Type | Required? | Description |
|-|-|-|-|
| Transcript annotations (combined_annotations.gtf) | GTF | Yes | Generated in Step 3, annotations of ORFs, peptides and transcripts |
| Transcript counts (bambu_transcript_counts.txt) | TXT/CSV | No | Generated in Step 1, transcript counts per sample |
| Peptide intensities | TXT | No | Generated in Step 2 (external), Peptide intensity data 'report.pr_matrix.tsv' |

**Note:** This tool is interactive and there is an option to download plots as PNG, JPEG, SVG and PDF.