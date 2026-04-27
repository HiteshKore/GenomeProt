# This script removes the redundant ORFs and considers longest ORF if ORF is part of longer ORF.
# It annotates them based on their location on the genome/transcript
# cd-hit (First install cd-hit commond line tool in linux (command to install:'conda install bioconda/label/cf201901::cd-hit'). Next, install py-cdhit library using  'pip install py-cdhit' command)
# Usage:python3 annotate_proteome.py gencode.vM33.chr_patch_hapl_scaff.annotation_chrX.gtf openprot_uniprotDb_mm.txt  ORFome_aa.txt proteome_database_transcripts.gtf <outdir> <canonical/all> <orf_length> <variant_protein_db/None> <organism:HUMAN,CAEEL,MOUSE,RAT,DROME,DANRE>
#######################################################################

import sys
import os
import re
from parse_reference_gtf import *
from annotate_proteome_functions import *

# function to remove intermediate files
def remove_file_if_exists(filename):
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except:
            pass

# function to identify amino acid residue changes in wildtype and variant protein sequence
def variant_protein_annotations(var_ORF_anno, var_orfs, wt_protein):
    res = ""
    for var_protein in var_orfs:
        similarity, aa_change = calculate_similarity(var_protein, wt_protein)   # pairwise alignment to calculate sequence similarity
        similarity = round(similarity, 2)
        if 96 < similarity < 100:
            sq_properties = calculate_sequence_properties(var_protein)
            coordinates, var_type = var_ORF_anno[var_protein].split('|')[:2]
            res = '\t'.join(map(str, [var_protein, aa_change, coordinates, var_type, sq_properties]))
    return res

# function to write ORF sequences into the proteome database FASTA file
def writeORFsIntoFASTA(orf_info_map, organism_info, proteomedb_filename):
    (organism, organism_latin_name) = organism_info
    fa_seqs = ''
    with open(proteomedb_filename, 'a') as f:
        for prot_sq, annotations in orf_info_map.items():
            accession = ""
            description = set()
            orf_coord = set()   # ORF coordinates
            GA = set()          # gene accession
            GN = set()          # gene name
            TA = set()          # transcript

            # protein_accession+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript+"|"+protein_description
            for header_info in annotations:
                cols = header_info.split('|')
                [accession, orf_coordinates, gene_id, gene_name, transcript] = cols[:5]
                orf_coord.add(orf_coordinates)
                GA.add(gene_id)
                GN.add(gene_name)
                TA.add(transcript)
                if len(cols) == 6:
                    description.add(cols[-1])

            orf_coord_s = ','.join(orf_coord)
            TA_s = ','.join(TA)
            GA_s = ','.join(GA)

            des_s = ','.join(description)
            if not des_s:
                des_s = '-'

            GN_s = ','.join(GN)
            if not GN_s:
                GN_s = '-'

            orf_type = "nv" if (len(description) == 0 and "ORF_" in accession) else "kn"

            fa_seqs += f">{orf_type}|{accession}|{accession}_{organism} {des_s} OS={organism_latin_name} CO={orf_coord_s} GA={GA_s} GN={GN_s} TA={TA_s}\n" + f"{prot_sq}\n"
        f.write(fa_seqs)

# function to annotate protein sequence based on genomic coordinates, biotypes, physioco-chemical properties
def get_protein_annotation(transcript, gene_id, gene_name, protein_des, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, orf_type, localisation, openprot_annotations, longest_orf, protein_status, orf_metadata_map, reading_frame_info, var_ORF_anno):
    # get reading frame information
    matching_indices = next((i for i, protein_seq_rf in enumerate(reading_frame_info) if protein_seq in protein_seq_rf), None)

    rf = reading_frame_info[matching_indices].split("|")[1]

    if transcript in var_transcript_ORF_map:
        var_prot_attr = variant_protein_annotations(var_ORF_anno, var_transcript_ORF_map[transcript], protein_seq)  # identify variants and calculate physico-chemical properties
        if var_prot_attr != "":
            var_seq = var_prot_attr.split("\t")[0]
            variants = var_prot_attr.split("\t")[1]
            var_seq_coordinates = var_prot_attr.split("\t")[2]
            chr = var_seq_coordinates.strip().split(":")[0]
            if strand == "-":
                orf_start = int(var_seq_coordinates.strip().split(":")[1].split("-")[0]) + 3
                orf_end = var_seq_coordinates.strip().split(":")[1].split("-")[1]
                rf = (int(orf_end) - 1) % 3
                var_seq_coordinates = chr + ":" + str(orf_start) + "-" + orf_end
            elif strand == "+":
                orf_start = var_seq_coordinates.strip().split(":")[1].split("-")[0]
                orf_end = int(var_seq_coordinates.strip().split(":")[1].split("-")[1]) - 3  # 3 nucleiotides substracted as ORFik counts stop codon position
                rf = (int(orf_start) - 1) % 3
                var_seq_coordinates = chr + ":" + orf_start + "-" + str(orf_end)

            var_type = var_prot_attr.split("\t")[3]  # variants
            var_seq_prop = var_prot_attr.split("\t")[4]  # sequence properties
            if var_type == "HM":  # homozygous variants
                orf_annotation = protein_accession + "_var" + "\t" + gene_id + "\t" + gene_name + "\t" + protein_des + "\t" + transcript + "\t" + strand + "\t" + transcript_biotype + "\t" + transcript_coordinates + "\t" + var_seq_coordinates + "\t" + str(rf) + "\t" + orf_type + "\t" + localisation + "\t" + openprot_annotations + "\t" + var_seq + "\t" + longest_orf + "\t" + protein_status + "\t" + variants + "\t" + calculate_sequence_properties(var_seq) + "\n"
                orf_metadata_map.setdefault(var_seq, []).append(orf_annotation)
            elif var_type == "HT":  # heterozygous variants
                orf_annotation = '\t'.join([protein_accession, gene_id, gene_name, protein_des, transcript, strand, transcript_biotype, transcript_coordinates, orf_coordinate, str(rf), orf_type, localisation, openprot_annotations, protein_seq, longest_orf, protein_status, '-', calculate_sequence_properties(protein_seq)]) + '\n'
                orf_metadata_map.setdefault(protein_seq, []).append(orf_annotation)  # wildtype sequence
                orf_annotation = protein_accession + "_var" + "\t" + gene_id + "\t" + gene_name + "\t" + protein_des + "\t" + transcript + "\t" + strand + "\t" + transcript_biotype + "\t" + transcript_coordinates + "\t" + var_seq_coordinates + "\t" + str(rf) + "\t" + orf_type + "\t" + localisation + "\t" + openprot_annotations + "\t" + var_seq + "\t" + longest_orf + "\t" + protein_status + "\t" + variants + "\t" + calculate_sequence_properties(var_seq) + "\n"
                orf_metadata_map.setdefault(var_seq, []).append(orf_annotation)  # variant sequence
    else:  # if transcript not in var_transcript_ORF_map
        orf_annotation = '\t'.join([protein_accession, gene_id, gene_name, protein_des, transcript, strand, transcript_biotype, transcript_coordinates, orf_coordinate, str(rf), orf_type, localisation, openprot_annotations, protein_seq, longest_orf, protein_status, '-', calculate_sequence_properties(protein_seq)]) + '\n'
        orf_metadata_map.setdefault(protein_seq, []).append(orf_annotation)

def main():
    args = sys.argv
    if 10 <= len(args) <= 12:
        print("Usage: python annotate_proteome.py <reference_gtf> <custom_openprot+uniprot_db> <ORFome_aa.txt> <ORFome_transcripts.gtf> <outdir> <canonical/all> <orf_length> <variant_protein_db/None> <organism> (num_threads) (memory_limit)")
        sys.exit(1)

    # Store the arguments
    [arg_reference_gtf_filename, arg_combined_protein_db_filename, arg_orfome_filename, arg_orfome_transcript_gtf_filename, arg_outdir, arg_canonical_or_all, arg_orf_length, arg_mutant_protein_db_filename, arg_organism] = args[1:10]

    # Check the arguments
    error_message = ""
    if not os.path.isfile(arg_reference_gtf_filename):
        error_message = f"The reference GTF file '{arg_reference_gtf_filename}' either does not exist or is not a file."
    elif not os.path.isfile(arg_combined_protein_db_filename):
        error_message = f"The combined protein database file '{arg_combined_protein_db_filename}' either does not exist or is not a file."
    elif not os.path.isfile(arg_orfome_filename):
        error_message = f"The ORFome annotation file '{arg_orfome_filename}' either does not exist or is not a file."
    elif not os.path.isfile(arg_orfome_transcript_gtf_filename):
        error_message = f"The ORFome transcripts GTF file '{arg_orfome_transcript_gtf_filename}' either does not exist or is not a file."
    elif os.path.exists(arg_outdir) and not os.path.isdir(arg_outdir):
        error_message = f"'{arg_outdir}' exists but it is not a directory."
    elif arg_mutant_protein_db_filename.lower() != "none" and not os.path.isfile(arg_mutant_protein_db_filename):
        error_message = f"The variant ORFome annotation file '{arg_mutant_protein_db_filename}' either does not exist or is not a file."
    elif arg_canonical_or_all.lower().strip() not in ["canonical", "all"]:
        error_message = f"The type of ORFs to be annotated can only be 'canonical' or 'all', not '{arg_canonical_or_all}' (case-insensitive)."

    if error_message:
        print(error_message)
        sys.exit(1)

    # If the output directory does not exist, try to create it
    if not os.path.exists(arg_outdir):
        try:
            os.mkdir(arg_outdir)
        except Exception as e:
            print(f"Failed to create the output directory '{arg_outdir}'. Error message: '{e}'.")
            sys.exit(1)

    # Parse the ORF length as a non-negative integer
    try:
        arg_orf_length = int(arg_orf_length)
        if arg_orf_length < 0:
            print(f"The ORF length must be provided as a non-negative integer.")
            sys.exit(1)
    except:
        print(f"The ORF length '{arg_orf_length}' cannot be parsed as an integer.")
        sys.exit(1)

    # Check if the organism is supported
    org_map = {"HUMAN": "Homo sapiens", "CAEEL": "Caenorhabditis elegans", "MOUSE": "Mus musculus", "RAT": "Rattus norvegicus", "DROME": "Drosophila melanogaster", "DANRE": "Danio rerio"}
    if arg_organism.upper().strip() not in org_map:
        print(f"The organism name '{arg_organism}' must be 'HUMAN', 'CAEEL', 'MOUSE', 'RAT', 'DROME' or 'DANRE' (case-insensitive).")
        sys.exit(1)
    organism_latin_name = org_map[arg_organism.upper().strip()]
    organism_info = (arg_organism.upper().strip(), organism_latin_name)

    proteomedb_filename = os.path.join(arg_outdir, "proteome_database.fasta")
    annotate_canonical_orfs_only = (arg_canonical_or_all.lower().strip() == "canonical")

    # Number of threads for cd-hit to use
    num_threads = 1
    try:
        num_threads = max(0, int(args[10]))
    except:
        pass

    # Maximum amount of memory (in MB) for cd-hit to use
    memory_limit = 800
    try:
        memory_limit = max(0, int(args[11]))
    except:
        pass

    # Load the combined OpenProt + UniProt annotation database
    openprot = {}   # OpenProt annotations: k:seq v:protein_id
    refprot = {}    # reference proteins annotated in uniprot/refseq/ensembl. k:seq v:protein id
    uniprot = {}    # UniProt annotations: k:seq v: trEMBL/reviewed|protein_id| gene description
    getDbAnnotations(arg_combined_protein_db_filename, openprot, refprot, uniprot)

    # Protein-coding genes
    protein_coding_gene_coordinates = {}    # k:chrom,v:set((gene start, gene end))

    # Transcripts
    transcript_biotypes = {}                # k:transcript_id i.e ENST,v: transcript biotype
    transcript_genome_coordinates = {}      # k:transcript_id i.e ENST,v: genome coordinates
    transcript_gene_id_map = {}             # k:transcript_id i.e ENST,v: gene id i.e. ENSG
    transcript_gene_name_map = {}           # k:transcript_id i.e ENST,v: gene name
    transcript_strand = {}                  # k:transcript_id i.e ENST,v: strand

    # UTRs and CDSs
    utr_coordinates = {}  # genome coordinates of transcript utrs k:transcript_id i.e ENST,k: genomic coordinates
    cds_coordinates = {}  # genome coordinates of transcript cds k:transcript_id i.e ENST,k: genomic coordinates

    with open(arg_reference_gtf_filename, 'r') as f:
        for raw_line in f:
            GP = GTFParser(raw_line)
            if not GP.valid:
                continue

            feature, strand, chrom, start, end, coordinates = GP.feature, GP.strand, GP.chrom, GP.start, GP.end, GP.coordinates
            gene_id, gene_name, gene_type = GP.gene_id, GP.gene_name, GP.gene_type
            transcript_id, transcript_type = GP.transcript_id, GP.transcript_type

            if feature == "gene" and gene_type == "protein_coding":     # Fetch coordinates for protein-coding genes
                protein_coding_gene_coordinates.setdefault(chrom, set()).add((start, end))
            elif feature == "transcript":
                transcript_biotypes[transcript_id] = transcript_type
                transcript_genome_coordinates[transcript_id] = coordinates
                transcript_strand[transcript_id] = strand
                transcript_gene_id_map[transcript_id] = gene_id
                transcript_gene_name_map[transcript_id] = gene_name

            if GP.featureExists("UTR"):
                utr_coordinates.setdefault(transcript_id, []).append(coordinates)
            if GP.featureExists("CDS"):
                cds_coordinates.setdefault(transcript_id, []).append(coordinates)

    with open(arg_orfome_transcript_gtf_filename, 'r') as f:
        for raw_line in f:
            if not ("BambuGene" in raw_line or "BambuTx" in raw_line):  # denovoGene denovoTx
                continue

            GP = GTFParser(raw_line)
            if not GP.valid:
                continue

            strand, coordinates = GP.strand, GP.coordinates
            gene_id, gene_name = GP.gene_id, GP.gene_name
            transcript_id = GP.transcript_id

            if GP.featureExists("transcript"):
                if "BambuTx" in transcript_id:
                    transcript_biotypes[transcript_id] = "novel"
                transcript_genome_coordinates[transcript_id] = coordinates
                transcript_strand[transcript_id] = strand
                transcript_gene_id_map[transcript_id] = gene_id
                transcript_gene_name_map[transcript_id] = gene_name

    var_transcript_ORF_map = {}
    var_ORF_anno = {}
    if arg_mutant_protein_db_filename.lower() != "none":    # variant file is optional
        with open(arg_mutant_protein_db_filename, 'r') as f:
            for raw_line in f:
                line = raw_line.strip()
                if line.startswith("transcript"):
                    continue
                columns = line.split('\t')[:4]
                transcript = columns[0]
                if transcript in transcript_gene_id_map and len(columns) == 4:  # the gencode version of reference GTF must be same
                    [var_sq, orf_coordinate, mutation_type] = columns[1:]
                    var_transcript_ORF_map.setdefault(transcript, []).append(var_sq)
                    var_ORF_anno[var_sq] = f"{orf_coordinate}|{mutation_type}"

    annotated_proteins = {}                 # store annotated proteins: k:protein sequence, v:transcript|orf_id|genomic_coordinates
    unannotated_proteins = {}               # stores unannotated proteins: k:protein sequence, v:orf_id|genomic_coordinates
    unannotated_protein_coordinates = {}    # stores unannotated protein coordinates as key: k: orf_id|genomic_coordinates,v:protein sequence
    transcript_orf_map = {}                 # stores all orfs for each transcript: k;transcript_id,v:protein_seq|orf_id
    reading_frame_map = {}

    # Load the ORFome annotation file
    with open(arg_orfome_filename, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith("ORF_id"):
                continue

            cols = [col.strip() for col in line.split('\t')]
            [orf_id, protein_seq, chrom, orf_start, orf_end, strand, reading_frame] = cols[:7]
            if (strand != '+') and (strand != '-'):
                continue

            transcript = orf_id.split('_')[0]
            reading_frame_map.setdefault(transcript, []).append(f"{protein_seq}|{reading_frame}")
            orf_start, orf_end = int(orf_start), int(orf_end)

            is_forward_strand = (strand == '+')
            orf_start += (3 if not is_forward_strand else 0)
            orf_end -= (3 if is_forward_strand else 0)          # 3 nucleotides subtracted as ORFik counts stop codon position
            orf_coordinate = f"{chrom}:{orf_start}-{orf_end}"

            if transcript in cds_coordinates:
                # check if transcript contains annotated CDS
                cds_start = int(cds_coordinates[transcript][0 if is_forward_strand else -1].split('-')[0].split(':')[1])
                cds_end = int(cds_coordinates[transcript][-1 if is_forward_strand else 0].split('-')[1])
                if orf_start >= cds_start and orf_end <= cds_end:
                    transcript_orf_map.setdefault(transcript, []).append(f"{protein_seq}|{orf_id}") # store all orfs of transcript
            else:
                transcript_orf_map.setdefault(transcript, []).append(f"{protein_seq}|{orf_id}")     # store all orfs of transcript

            # separate annnotated and unannotated proteins
            if protein_seq in uniprot:
                uniprot_accession = uniprot[protein_seq].split('|')[1]
                annotated_proteins.setdefault(protein_seq, []).append(f"{uniprot_accession}|{orf_id}|{orf_coordinate}")
            elif protein_seq in refprot:
                uniprot_ids = uniprot.values()
                refprot_id = refprot[protein_seq]
                is_id_found = any(refprot_id in uniprot_id for uniprot_id in uniprot_ids)
                if not is_id_found:  # avoid incorrect protein ids assigned in openprot
                    annotated_proteins.setdefault(protein_seq, []).append(f"{refprot_id}|{orf_id}|{orf_coordinate}")
            else:
                unannotated_proteins.setdefault(protein_seq, []).append(f"{orf_id}|{orf_coordinate}")
                unannotated_protein_coordinates[f"{orf_id}|{orf_coordinate}"] = protein_seq

    # Find the longest ORF in the transcripts
    transcript_longest_orf_map = {}  # stores longest ORF in transcript: k: orf_id v:protein_seq
    for transcript_orf_recs in transcript_orf_map.values():
        longest_orf_seq, longest_orf_id, longest_orf_len = '', '', -1
        for orf_rec in transcript_orf_recs:
            orf_seq, orf_id = [col.strip() for col in orf_rec.split('|')[:2]]
            orf_len = len(orf_seq)
            if orf_len > longest_orf_len:
                longest_orf_seq, longest_orf_id, longest_orf_len = orf_seq, orf_id, orf_len
        transcript_longest_orf_map[longest_orf_id] = longest_orf_seq

    # Store protein sequences in a file for cd-hit clustering
    orf_temp_file = os.path.join(arg_outdir, "orf_temp.txt")
    with open(orf_temp_file, 'w') as f:
        entries = ''
        for protein_seq, protein_info in unannotated_proteins.items():
            entries += f">{protein_info[0]}\n" + f"{protein_seq}\n"
        for protein_seq, protein_info in annotated_proteins.items():
            entries += f">{protein_info[0]}\n" + f"{protein_seq}\n"
        _ = f.write(entries)

    # Perform cd-hit clustering on unannotated proteins to consider longest representative protein sequence
    orf_clus_file = os.path.join(arg_outdir, "cdhit_out")
    SeqClust(orf_temp_file, orf_clus_file, num_threads, memory_limit)

    # annotations of known proteins
    orf_metadata_annotated_proteins = {}  # key:protein_seq, val:metadata
    metadata_records = []

    for protein_seq, annotated_protein_infos in annotated_proteins.items():
        is_protein_seq_in_uniprot = (protein_seq in uniprot)
        for annotated_protein_info in annotated_protein_infos:
            [orf_id, orf_coordinate] = [col.strip() for col in annotated_protein_info.split('|')[1:3]]
            transcript = orf_id.split('_')[0]
            if transcript not in transcript_biotypes:
                continue

            gene_name = transcript_gene_name_map[transcript]
            if not gene_name:
                gene_name = '-'

            gene_id = transcript_gene_id_map[transcript]
            strand = transcript_strand[transcript]
            transcript_biotype = transcript_biotypes[transcript]
            transcript_coordinates = transcript_genome_coordinates[transcript]
            longest_orf = 'Y' if orf_id in transcript_longest_orf_map else 'N'

            protein_status = '-'
            protein_description = '-'
            if is_protein_seq_in_uniprot:   # UniProt proteins
                [database, protein_accession, protein_description] = uniprot[protein_seq].split('|')[:3]
                protein_status = "reviewed(Swiss-Prot)" if database == "sp" else "unreviewed(TrEMBL)"
            else:                           # RefProt proteins
                protein_accession = refprot[protein_seq]

            # annotate ORFs
            get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "annotated", "CDS", "-", longest_orf, protein_status, orf_metadata_annotated_proteins, reading_frame_map[transcript], var_ORF_anno)

    # Calculate annotated ORF outputs
    annotated_orf_map = {}  # k:protein_seq, v:fasta header annotations
    for protein_seq, annotations in orf_metadata_annotated_proteins.items():
        for protein_anno in annotations:
            cols = protein_anno.split('\t')
            [protein_accession, gene_id, gene_name, protein_description, transcript] = cols[:5]
            orf_coordinate, longest_orf = cols[8], cols[14]

            if annotate_canonical_orfs_only and longest_orf != 'Y':
                continue

            metadata_records.append(protein_anno)
            fa_header = '|'.join([protein_accession, orf_coordinate, gene_id, gene_name, transcript, protein_description])
            annotated_orf_map.setdefault(protein_seq, []).append(fa_header)

    # Write annotated ORFs into the proteome database FASTA file
    writeORFsIntoFASTA(annotated_orf_map, organism_info, proteomedb_filename)

    # Novel proteins
    orf_annotation_map = {}  # stores ORF annotations k:temporory ORF_id, annotations
    cdhit_clustering_output_file = os.path.join(arg_outdir, "cdhit_out.clstr")  # cluster file
    with open(cdhit_clustering_output_file, 'r') as cdhit_clustering_output:
        counter = 1     # Counter for temporary ORF IDs
        protein_description = "-"
        protein_status = "-"
        for raw_line in cdhit_clustering_output:
            line = raw_line.strip()
            if not ('*' in line and (">Bambu" in line or ">ENS" in line)):
                continue

            # remove representative sequences of known proteins
            if re.search(r">ENS\w*P\d{11}", line):
                continue

            longest_seq_orf_id = line.split('>')[1].split("...")[0]                 # ORF ID of the longest representative ORF in the cluster
            protein_seq = unannotated_protein_coordinates[longest_seq_orf_id]       # Protein sequence
            openprot_id = openprot[protein_seq] if protein_seq in openprot else '-' # Check if the protein is annotated in OpenProt

            protein_accession = f"ORF_{counter}"
            counter += 1

            for unannotated_protein_info in unannotated_proteins[protein_seq]:  # accessing all ORF coordinates for a given protein sequence
                [orf_id, orf_coordinate] = unannotated_protein_info.split('|')[:2]
                transcript = orf_id.split('_')[0]
                if transcript not in transcript_biotypes:
                    continue

                transcript_biotype = transcript_biotypes[transcript]
                if transcript_biotype == "protein_coding" and transcript not in utr_coordinates:
                    continue

                transcript_coordinates = transcript_genome_coordinates[transcript]
                gene_id = transcript_gene_id_map[transcript]
                strand = transcript_strand[transcript]
                longest_orf = 'Y' if orf_id in transcript_longest_orf_map else 'N'

                gene_name = transcript_gene_name_map[transcript]
                if not gene_name:
                    gene_name = '-'

                if transcript_biotype == "protein_coding":  # protein-coding transcripts
                    localisation = UTRAnnotations(utr_coordinates[transcript], orf_coordinate, strand, cds_coordinates[transcript])
                    if not ("UTR" in localisation and localisation != "CDS:3UTR"):  # ORF overlaps with UTR region
                        continue
                    if len(protein_seq) < arg_orf_length:   # for short uORFs, remove CDS annotation
                        localisation = localisation.split(':')[0]
                elif transcript in utr_coordinates:         # transcripts with UTRs but are 'non-coding'
                    localisation = UTRAnnotations(utr_coordinates[transcript], orf_coordinate, strand, cds_coordinates[transcript])
                    if "UTR" in localisation and localisation != "CDS:3UTR":        # ORF overlaps with UTR region
                        if len(protein_seq) < arg_orf_length:
                            localisation = localisation.split(':')[0]
                    else:   # if ORF is not in UTR regions
                        is_gene_overlap = isIntergenic(orf_coordinate, protein_coding_gene_coordinates)
                        localisation = "gene_overlap" if is_gene_overlap else "intergenic"
                else:                                       # non-coding transcripts that don't have UTRs
                    is_gene_overlap = isIntergenic(orf_coordinate, protein_coding_gene_coordinates)
                    localisation = "gene_overlap" if is_gene_overlap else "intergenic"

                # annotate ORFs
                get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", localisation, openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)

    # Compute the output of novel ORFs
    novel_orf_map = {}          # Replace the temporary ORF IDs based on their index
    novel_orf_id_dict = {}      # Store the IDs of novel ORFs
    novel_orf_id_counter = 1
    variant_suffix = "_var"

    for protein_seq, annotations in orf_annotation_map.items():
        for protein_annotation in annotations:
            cols = protein_annotation.split('\t')
            [accession, gene_id, gene_name], transcript, orf_coordinate, longest_orf = cols[:3], cols[4], cols[8], cols[14]
            if annotate_canonical_orfs_only and longest_orf != 'Y':
                continue

            temp_accession = accession.replace(variant_suffix, '')
            novel_orf_num = novel_orf_id_dict.get(temp_accession, novel_orf_id_counter)
            if novel_orf_num == novel_orf_id_counter:
                novel_orf_id_dict[temp_accession] = novel_orf_id_counter
                novel_orf_id_counter += 1

            new_accession = f"ORF_{novel_orf_num}" + (variant_suffix if variant_suffix in accession else '')
            fa_header = '|'.join([new_accession, orf_coordinate, gene_id, gene_name, transcript])
            novel_orf_map.setdefault(protein_seq, []).append(fa_header)
            revised_annotation = protein_annotation.replace(accession, new_accession)
            metadata_records.append(revised_annotation)

    metadata_records_uniq = list(dict.fromkeys(metadata_records))

    # Write the proteome database metadata file
    proteome_database_metadata_filename = os.path.join(arg_outdir, "proteome_database_metadata.txt")
    with open(proteome_database_metadata_filename, 'w') as f:
        f.write('\t'.join(["accession", "gene", "gene_symbol", "protein_description", "transcript", "strand", "transcript_biotype", "transcript_coordinates", "orf_genomic_coordinates", "reading_frame", "orf_type", "localisation", "openprot_id", "protein_sequence", "longest_orf_in_transcript", "uniprot_status", "amino_acid_change", "molecular_weight(kDA)", "isoelectric_point", "hydrophobicity", "aliphatic_index"]) + '\n')
        f.writelines(metadata_records_uniq)

    # Write novel ORFs into the proteome database FASTA file
    writeORFsIntoFASTA(novel_orf_map, organism_info, proteomedb_filename)

    # Remove intermediate files
    remove_file_if_exists(orf_temp_file)
    remove_file_if_exists(orf_clus_file)
    remove_file_if_exists(cdhit_clustering_output_file)

if __name__ == "__main__":
    main()