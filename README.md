# GenomeProt: an integrated proteogenomics analysis platform for long-read RNA-Seq datasets


## Contents

- [Installation](#installation)
- [General usage](#general-usage)
  - [Database generation](#1-database-generation)
  - [Proteomics](#2-proteomics)
  - [Integration](#3-integration)
  - [Visualisation](#4-visualisation)
- [Detailed input and output descriptions](#detailed-input-and-output-descriptions)
  - [Database generation](#1-database-generation-1)
  - [Proteomics](#2-proteomics-1)
  - [Integration](#3-integration-1)
  - [Visualisation](#4-visualisation-1)

## Installation

### Option 1 (recommended): Access GenomeProt online
https://genomeprot.researchsoftware.unimelb.edu.au/

### Option 2: Run the shiny application with Docker
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```
docker run --rm -p 3838:3838 josieg/genomeprot:v1
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it’s stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.

### Option 3: Locally install the shiny application

The application has substantial dependencies that we have provided as a conda environment file. 

Clone this repository:
```
git clone https://github.com/josiegleeson/GenomeProt.git
```

Build the conda environment:
```
cd GenomeProt
Rscript install.R conda_env.yaml
```

Unzip the uniprot+openprot reference file in the GenomeProt/data directory.
```
cd GenomeProt/GenomeProt/data
unzip openprot_uniprotDb_hs.txt.zip
unzip *
```

Activate the environment and then run the app from the command line:
```
conda activate GenomeProt_env
Rscript -e "shiny::runApp('path/to/app/GenomeProt/', host='0.0.0.0', port=3838)"
```


## General usage 

GenomeProt is an integrated proteogenomics platform with four modules: 1) database generation, 2) proteomics (under development), 3) integration, and 4) visualisation.

### 1. Database generation

The first module generates a custom proteome database to perform proteomics searches. The module accepts RNA sequencing FASTQ files, BAM files or GTF annotation files from both short-read and long-read sequencing platforms (Illumina, Oxford Nanopore and PacBio). The main output from this module is a FASTA file with candidate protein sequences and a metadata file with details of each candidate protein. 

For long-read data, transcript discovery is performed with Bambu. For short-read data, no discovery steps are performed, transcripts are instead directly quantified based on the reference transcriptome using Salmon. GenomeProt currently supports open reading frame (ORF)  identification and database generation for six model organisms: fruit fly, roundworm, zebrafish, rat, mouse, and human. Users can specify an option to include short upstream ORF (uORF) and downstream ORF (dORF) protein sequences >10 amino acids (AA). Protein sequences are generated based on a user defined minimum length set to >30 AA by default. Users can also optionally provide a VCF file to incorporate single nucleotide variants (SNVs) into the genome to generate variant protein sequences. 

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

### 2. Proteomics

Currently under development. Please use FragPipe to process proteomics data with the custom database.

### 3. Integration

This module integrates proteomics and transcriptomics data. Peptides are associated back to transcript isoforms and mapped to spliced genomic coordinates for downstream visualisation. This generates BED12 file of transcripts, ORFs  and peptides for visualisation in the UCSC genome browser and a combined GTF file for visualisation within the app. An html report is also created that provides a summary of identified known and novel transcripts, uniquely mapping peptides, and known and novel ORFs.

#### Key outputs:

- Peptide-to-transcript mappings with spliced genomic coordinates. 
- BED12 files for visualisation in UCSC genome browser. 
- HTML report summarising identified transcripts, peptides and ORFs. 

### 4. Visualisation

There are currently two visualisation options included. The first (4a) is a custom tool using ggtranscript that generates peptide mapping plots along transcript isoforms with quantitative peptide intensities and transcript expression data. This allows users to visualise transcript and peptide abundance across different experimental conditions. This module requires the combined GTF file generated in the integration step, and optionally inputs transcript counts from Module 1 and peptide intensities from external proteomics analysis. The advantage of this tool is the included gene filtering options to quickly search for features of interest. The other tool (4b) included is IsoVis, a webserver that can be accessed within the app or externally. This is a more comprehensive tool that requires the `transcripts_and_ORFs_for_isovis.gtf` file and optionally the peptide intensities and transcript expression data.

#### Features:

- Requires `combined_annotations.gtf` or `transcripts_and_ORFs_for_isovis.gtf` (from Module 3).
- Optionally input transcript isoform counts and peptide intensities.
- Allows export of plots as PDFs. 

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


| Output  | File   | File Type | Description      |
|-------------------------------------|------------------------------------|------------|---------------------------------------------------------------------------------------------|
| Database   | proteome_database.fasta  | FASTA | Amino acid sequences of all ORFs in the data  |
| Database metadata | proteome_database_metadata.txt | TXT  | Information on each ORF in the data    |
| Database transcripts (long-read data only) | proteome_database_transcripts.gtf | GTF  | Annotations of transcripts used to generate the database   |
| Transcript counts (FASTQ or BAM input only) | transcript_counts.txt  | TXT  | Transcript count file with samples as columns and transcripts as rows  |
| Bambu transcript class codes (long-read FASTQ or BAM input only) | novel_transcript_classes.csv  | CSV  | Bambu transcript classification (see Bambu documentation)    |
| GFFcompare transcript class codes (long-read FASTQ or BAM input only) | gffcompare.tmap.txt | TXT  | GFFcompare transcript classification (see GFFcompare documentation) |
| Genome aligned reads (long-read FASTQ only) | sample1.bam, sample2.bam   | BAM  | Genome aligned reads. Only output if the input was FASTQ reads    |


Proteome FASTA examples and header formats: 
```
>protein_accession|CO=genomic_coordinates GA=gene_accession GN=gene_name TA=transcript_accession
MCGNNMSAPMPAVVPAARKATAAVIFLHGLGDTGHGWAEAFAGIKSPHIKYICPHAPVMPVTLNMNMAMPSWFDIVGLSPDSQEDESGIKQAAETVKALIDQEVKNGIPSNRIILGGFSQGPINSANRDISVLQCHGDCDPLVPLMFGSLTVERLKALINPANVTFKIYEGMMHSSCQQEMMDVKHFIDKLLPPID
>P10711|CO=chr1:4928137-4966584 GA=ENSMUSG00000033813.16 GN=Tcea1 TA=ENSMUST00000081551.14
MEDEVVRIAKKMDKMVQKKNAAGALDLLKELKNIPMTLELLQSTRIGMSVNALRKQSTDEEVTSLAKSLIKSWKKLLDGPSTDKDPEEKKKEPAISSQNSPEAREESSSSSNVSSRKDETNARDTYVSSFPRAPSTSDSVRLKCREMLAAALRTGDDYVAIGADEEELGSQIEEAIYQEIRNTDMKYKNRVRSRISNLKDAKNPNLRKNVLCGNIPPDLFARMTAEEMASDELKEMRKNLTKEAIREHQMAKTGGTQTDLFTCGKCKKKNCTYTQVQTRSADEPMTTFVVCNECGNRWKFC
>ORF_3|CO=chr2:53029193-53081430 GA=ENSMUSG00000061136.17 GN=Prpf40a TA=ENSMUST00000209364.3
MQATPSEAGGESPQSCLSVSRSDWTVGKPVSLLAPLIPPRSSGQPLPFGPGGRQPLRSLLVGMCSGSGRRRSSLSPTMRPGTGAERGGLMMGHPGMHYAPMGMHPMGQRANMPPVPHGMMPQMMPPMGG
```
**Note:** Unannotated ORFs are denoted by "ORF_" and variant proteins by “mORF_” followed by a unique number. UniProt or RefSeq accessions are retained for annotated proteins.

#### Open reading frame (ORF) category definitions:

| Type  | Definition    |
|----------------|-----------------------------------------------------------------------------------------------------|
| CDS  | Annotated in UniProt or RefSeq    |
| 5UTR  | Coordinates are within the 5’ UTR region of an mRNA transcript    |
| 3UTR  | Coordinates are within the 3’ UTR region of an mRNA transcript    |
| 5UTR:CDS  | Start site is within the 5’ UTR region and stop site is within the CDS region of an mRNA transcript |
| gene_overlap  | Encoded by a transcript that overlaps a region with annotated protein-coding genes   |
| intergenic | Encoded by a transcript that does not overlap a region with annotated protein-coding genes   |


### 2. Proteomics

We recommend installing and running [FragPipe](https://github.com/Nesvilab/FragPipe) for analysing mass spectrometry-based proteomics data.

| Input  | File Type | Required? | Description     |
|------------------------------------|-----------|-----------|-----------------------------------------------------------------------------------------------|
| Mass spec data  | mzML, RAW  | Yes  | Mass spec files  |
| Database (proteome_database.fasta) | FASTA | Yes  | Generated in Module 1. Amino acid sequences of all ORFs in the data  |

The output file generated is typically `peptides.txt` or `report.pr_matrix.tsv`.


### 3. Integration

| Input  | File Type | Required? | Description     |
|------------------------------------|-----------|-----------|-----------------------------------------------------------------------------------------------|
| Proteomics peptide data  | TSV/TXT  | Yes  | Peptide results. Typically, 'peptides.txt', ‘peptide.tsv’ or ‘report.pr_matrix.tsv’  |
| Database (proteome_database.fasta) | FASTA | Yes  | Generated in Module 1. Amino acid sequences of all ORFs in the data  |
| Database metadata (proteome_database_metadata.txt) | TXT  | Yes  | Generated in Module 1. Information on each ORF in the data   |
| Database transcripts (proteome_database_transcripts.gtf) | GTF  | Yes  | Generated in Module 1. Annotations of transcripts used to generate the database   |


| Output | File  | File Type | Description    |
|------------------------------|--------------------------|-----------|----------------------------------------------------------------|
| Peptide information   | peptide_info.csv   | CSV  | Main results file with peptide mapping data |
| Report | summary_report.html | HTML | Summary report  |
| Combined annotation data | combined_annotations.gtf | GTF  | Annotations of peptides, ORFs, and transcripts for visualisation |
| Transcript and ORF annotation data | transcripts_and_ORFs_for_isovis.gtf | GTF  | Annotations of ORFs and transcripts for visualisation with IsoVis |
| Peptide coordinates   | peptides.bed12  | BED12 | Peptide spliced genomic coordinates   |
| ORF coordinates   | ORFs.bed12 | BED12 | ORF spliced genomic coordinates   |
| Transcript coordinates  | transcripts.bed12  | BED12 | Transcript spliced genomic coordinates |


#### Description of `peptides_info.csv` output:

| Column Name | Description     | Class  |
|--------------------------------|----------------------------------------------------------------------------------------|--------------|
| peptide | Peptide sequence   | character  |
| accession  | Protein accession    | character  |
| PID   | protein_accession\|CO=genomic_coordinates (included for compatibility with FASTA header) | character  |
| transcript_id   | ENSEMBL or novel transcript ID     | character  |
| gene_id | ENSEMBL gene ID   | character  |
| gene_name  | Gene name/symbol   | character  |
| strand  | Strand (+ or -)   | character  |
| number_exons | Number of exons spanned by the peptide   | integer |
| transcript_length   | Transcript length (nt)     | integer |
| transcript_biotype   | Transcript biotype from Gencode    | character  |
| simplified_biotype   | Simplified transcript biotype    | character  |
| protein_length   | Protein length (AA)   | integer |
| orf_genomic_coordinates  | Genomic coordinates of ORF    | numeric |
| orf_type | Annotated, unannotated, or variant protein   | character  |
| localisation | ORF location in the genome (see ORF category definitions)  | character  |
| uniprot_status   | Review status in UniProt (reviewed/unreviewed)    | character  |
| openprot_id | OpenProt ID if present     | character  |
| molecular_weight(kDA)   | Molecular weight of protein (KDa)   | numeric |
| isoelectric_point   | Isoelectric point of ORF calculated using pKa scale EMBOSS (Rice et al., 2000)   | numeric |
| hydrophobicity   | Hydrophobicity profile of ORF calculated using Kyte-Doolittle scale (Kyte et al., 1982) | numeric |
| aliphatic_index | Aliphatic index of ORF (Ikai 1980)  | numeric |
| longest_orf_in_transcript | Longest ORF in the transcript (longest within CDS regions for known proteins)   | true/false  |
| peptide_ids_gene | Is peptide uniquely mapped to gene?    | true/false  |
| peptide_ids_orf | Is peptide uniquely mapped to ORF?    | true/false  |
| peptide_ids_transcript   | Is peptide uniquely mapped to transcript?   | true/false  |
| shared_novel_protein_peptide  | Is peptide shared with other novel proteins?  | true/false  |
| orf_identified   | Is ORF identified with unique peptide evidence?  | true/false  |
| gene_identified | Is gene identified with unique peptide evidence?  | true/false  |
| transcript_identified   | Is transcript identified with unique peptide evidence?    | true/false  |


### 4. Visualisation

#### 4a. Visualisation: ggtranscript

| Input  | File Type | Required? | Description   |
|--------------------------------|-----------|-----------|---------------------------------------------------------|
| Combined annotations (combined_annotations.gtf)  | GTF  | Yes  | Generated in Module 3, annotations of peptides, ORFs, and transcripts    |
| Transcript counts (transcript_counts.txt)   | TXT/CSV  | No  | Generated in Module 1, transcript counts per sample    |
| Peptide intensities  | TXT  | No  | Peptide intensity data ‘report.pr_matrix.tsv’  |

**Note:** There is an option to download plots as a PDF.

#### 4b. Visualisation: IsoVis

| Input  | File Type | Required? | Description   |
|--------------------------------|-----------|-----------|---------------------------------------------------------|
| Transcript annotations (transcripts_and_ORFs_for_isovis.gtf)  | GTF  | Yes  | Generated in Module 3, annotations of ORFs, and transcripts    |
| Peptide track (peptides.bed12)  | BED12  | Yes  | Generated in Module 3, coordinates of each peptide    |
| Transcript counts (transcript_counts.txt)   | TXT/CSV  | No  | Generated in Module 1, transcript counts per sample    |
| Peptide intensities  | TXT  | No  | Peptide intensity data ‘report.pr_matrix.tsv’  |

**Note:** This tool is interactive and there is an option to download plots as a PDF.
