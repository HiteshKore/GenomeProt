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

def main():
    args = sys.argv
    if len(args) != 10:
        print(f"Usage: {args[0]} annotate_proteome.py <reference_gtf> <custom_openprot+uniprot_db> <ORFome_aa.txt> <ORFome_transcripts.gtf> <outdir> <canonical/all> <orf_length> <variant_protein_db/None> <organism>")
        sys.exit(1)

    # Store the arguments
    [arg_reference_gtf_filename, arg_combined_protein_db_filename, arg_orfome_filename, arg_orfome_transcript_gtf_filename, arg_outdir, arg_canonical_or_all, arg_orf_length, arg_mutant_protein_db_filename, arg_organism] = args[1:]

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
    elif arg_mutant_protein_db_filename != "None" and not os.path.isfile(arg_mutant_protein_db_filename):
        error_message = f"The ORFome annotation file '{arg_mutant_protein_db_filename}' either does not exist or is not a file."

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

    fw_proteomedb = open(os.path.join(arg_outdir, "proteome_database.fasta"), 'w')

    # custom openprot + uniprot annotation database
    openprot = {}   # comprise OpenProt annotations: k:seq v:protein_id
    refprot = {}    # comprise reference proteins annotated in uniprot/refseq/ensembl. k:seq v:protein id
    uniprot = {}    # comprise UniProt annotations: k:seq v: trEMBL/reviewed|protein_id| gene description
    getDbAnnotations(arg_combined_protein_db_filename, openprot, refprot, uniprot)  # OpenProt annotations

    # Gene
    gene_biotype = {}  # k:gene_id i.e ENSG,k:gene biotype
    gene_coordinates = {}  # k:gene_id,k:gene coordinates
    protein_coding_gene_coordinates = {}  # k:gene_id,k:gene coordinates

    # transcripts
    transcript_biotypes = {}  # k:transcript_id i.e ENST,k: transcript biotype
    transcript_genome_coordinates = {}  # k:transcript_id i.e ENST,k: genome coordinates
    transcript_gene_id_map = {}  # k:transcript_id i.e ENST,k: gene id i.e. ENSG
    transcript_gene_name_map = {}  # k:transcript_id i.e ENST,k: gene name
    transcript_strand = {}  # k:transcript_id i.e ENST,k: strand
    # exon
    exon_lengths = {}  # transcript exon lengths k:transcript_id i.e ENST,k: lenghts of exons
    exon_coordinates = {}  # genome coordinates of transcript exons k:transcript_id i.e ENST,k: genomic coordinates

    # utr
    utr_coordinates = {}  # genome coordinates of transcript utrs k:transcript_id i.e ENST,k: genomic coordinates
    # cds coordinates
    cds_coordinates = {}  # genome coordinates of transcript cds k:transcript_id i.e ENST,k: genomic coordinates

    with open(arg_reference_gtf_filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            GP = GTFParser(line.strip())
            if GP.feature == "gene":
                gene_biotype[GP.ensgene] = GP.genetype
                gene_coordinates[GP.ensgene] = GP.getCoordinates()
                if GP.genetype == "protein_coding":  # Feching coordinates for protein coding genes
                    protein_coding_gene_coordinates[GP.ensgene] = GP.getCoordinates()
            if GP.feature == "transcript":
                transcript_biotypes[GP.transcriptid] = GP.transcript_type
                transcript_genome_coordinates[GP.transcriptid] = GP.getCoordinates()
                transcript_strand[GP.transcriptid] = GP.strand
                transcript_gene_id_map[GP.transcriptid] = GP.ensgene
                transcript_gene_name_map[GP.transcriptid] = GP.genename

            if GP.featureExists("exon"):
                exon_lengths.setdefault(GP.transcriptid, []).append(GP.end - GP.start + 1)  # stores exon length
                exon_coordinates.setdefault(GP.transcriptid, []).append(str(GP.start) + "-" + str(GP.end))  # stores exon gen
            if GP.featureExists("UTR"):
                utr_coordinates.setdefault(GP.transcriptid, []).append(GP.getCoordinates())
            if GP.featureExists("CDS"):
                cds_coordinates.setdefault(GP.transcriptid, []).append(GP.getCoordinates())

    # ORFome_transcripts gtf
    with open(arg_orfome_transcript_gtf_filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if "BambuGene" in line.strip() or "BambuTx" in line.strip():  # denovoGene denovoTx
                GP = GTFParser(line.strip())
                if GP.featureExists("exon"):
                    exon_lengths.setdefault(GP.transcriptid, []).append(GP.end - GP.start + 1)  # stores exon length
                    exon_coordinates.setdefault(GP.transcriptid, []).append(str(GP.start) + "-" + str(GP.end))  # stores exon gen
                if GP.featureExists("transcript"):
                    transcript_genome_coordinates[GP.transcriptid] = GP.getCoordinates()
                    transcript_strand[GP.transcriptid] = GP.strand
                    transcript_gene_id_map[GP.transcriptid] = GP.ensgene
                    transcript_gene_name_map[GP.transcriptid] = GP.genename
                    if "BambuTx" in GP.transcriptid:
                        transcript_biotypes[GP.transcriptid] = "novel"
                    if "BambuGene" in GP.genename:
                        gene_biotype[GP.genename] = "novel"

    var_transcript_ORF_map = {}
    var_ORF_anno = {}
    if arg_mutant_protein_db_filename != "None" and os.path.isfile(arg_mutant_protein_db_filename):  # variant file is optional
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

    # function to annotate protein sequence based on genomic coordinates, biotypes, physioco-chemical properties
    def get_protein_annotation(transcript, gene_id, gene_name, protein_des, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, orf_type, localisation, openprot_annotations, longest_orf, protein_status, orf_metadata_map, reading_frame_info, var_ORF_anno):
        # get reading frame information
        matching_indices = next((i for i, protein_seq_rf in enumerate(reading_frame_info) if protein_seq in protein_seq_rf), None)

        rf = reading_frame_info[matching_indices].split("|")[1]

        if transcript in var_transcript_ORF_map.keys():
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
                    orf_annotation = protein_accession + "\t" + gene_id + "\t" + gene_name + "\t" + protein_des + "\t" + transcript + "\t" + strand + "\t" + transcript_biotype + "\t" + transcript_coordinates + "\t" + orf_coordinate + "\t" + str(rf) + "\t" + orf_type + "\t" + localisation + "\t" + openprot_annotations + "\t" + protein_seq + "\t" + longest_orf + "\t" + protein_status + "\t-\t" + calculate_sequence_properties(protein_seq) + "\n"
                    orf_metadata_map.setdefault(protein_seq, []).append(orf_annotation)  # wildtype sequence
                    orf_annotation = protein_accession + "_var" + "\t" + gene_id + "\t" + gene_name + "\t" + protein_des + "\t" + transcript + "\t" + strand + "\t" + transcript_biotype + "\t" + transcript_coordinates + "\t" + var_seq_coordinates + "\t" + str(rf) + "\t" + orf_type + "\t" + localisation + "\t" + openprot_annotations + "\t" + var_seq + "\t" + longest_orf + "\t" + protein_status + "\t" + variants + "\t" + calculate_sequence_properties(var_seq) + "\n"
                    orf_metadata_map.setdefault(var_seq, []).append(orf_annotation)  # variant sequence

        else:  # if transcript not in var_transcript_ORF_map
            orf_annotation = protein_accession + "\t" + gene_id + "\t" + gene_name + "\t" + protein_des + "\t" + transcript + "\t" + strand + "\t" + transcript_biotype + "\t" + transcript_coordinates + "\t" + orf_coordinate + "\t" + str(rf) + "\t" + orf_type + "\t" + localisation + "\t" + openprot_annotations + "\t" + protein_seq + "\t" + longest_orf + "\t" + protein_status + "\t-\t" + calculate_sequence_properties(protein_seq) + "\n"
            orf_metadata_map.setdefault(protein_seq, []).append(orf_annotation)

    # function to generate fasta sequence for proteome database
    def formatFastaheader(orf_info_map, outfile):  # takes annotated ORFs map and output fasta file handler
        for prot_sq, annotations in orf_info_map.items():
            accession = ""
            org = arg_organism
            orgmap = {"HUMAN": "Homo sapiens", "CAEEL": "Caenorhabditis elegans", "MOUSE": "Mus musculus", "RAT": "Rattus norvegicus", "DROME": "Drosophila melanogaster", "DANRE": "Danio rerio"}

            orf_coord = []  # ORF coordinates
            GA = []  # gene accession
            GN = []  # gene name
            TA = []  # transcript\
            description = []
            # protein_accession+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript+"|"+protein_description
            for header_info in annotations:
                accession = header_info.split("|")[0]  # accession
                orf_coord.append(header_info.split("|")[1])  # ORF coordinates
                GA.append(header_info.split("|")[2])  # gene accession
                GN.append(header_info.split("|")[3])  # gene name
                TA.append(header_info.split("|")[4])  # transcript
                if len(header_info.split("|")) == 6:
                    description.append(header_info.split("|")[5])
            orf_coord_s = "CO=" + ",".join(list(set(orf_coord)))
            TA_s = "TA=" + ",".join(list(set(TA)))
            GA_s = "GA=" + ",".join(list(set(GA)))
            des_s = ",".join(list(set(description)))

            if len(GN) == 0:  # varify
                GN.append("-")
            GN_s = "GN=" + ",".join(list(set(GN)))
            if len(list(set(description))) > 0:
                fa_seq = ">kn|" + accession + "|" + accession + "_" + org + " " + des_s + " OS=" + orgmap[org] + " " + orf_coord_s + " " + GA_s + " " + GN_s + " " + TA_s + "\n" + prot_sq + "\n"
            else:
                if "ORF_" in accession:
                    fa_seq = ">nv|" + accession + "|" + accession + "_" + org + " -" + " OS=" + orgmap[org] + " " + orf_coord_s + " " + GA_s + " " + GN_s + " " + TA_s + "\n" + prot_sq + "\n"
                else:
                    fa_seq = ">kn|" + accession + "|" + accession + "_" + org + " -" + " OS=" + orgmap[org] + " " + orf_coord_s + " " + GA_s + " " + GN_s + " " + TA_s + "\n" + prot_sq + "\n"

            outfile.write(fa_seq)

    # ORFomedb
    unannotated_proteins = {}  # stores unannotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    unannotated_protein_coordinates = {}  # stores unannotated protein coordinates as key: k: transcript|orf_number|genomic_coordinates,v:protein sequence
    transcript_orf_map = {}  # stores all orfs for each transcript: k;transcript_id,v:protein_seq|orf_id

    annotated_proteins = {}  # store annotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    reading_frame_map = {}

    with open(arg_orfome_filename, 'r') as f:
        for i in f:
            if i.startswith("ORF_id"):
                continue
            orf_id = i.strip().split("\t")[0].strip()
            protein_seq = i.strip().split("\t")[1].strip()
            transcript = orf_id.split("_")[0]
            reading_frame = i.strip().split("\t")[6]
            reading_frame_map.setdefault(transcript, []).append(protein_seq + "|" + reading_frame)
            strand = i.strip().split("\t")[5].strip()
            chr = i.strip().split("\t")[2].strip()
            if strand == "-":
                orf_start = int(i.strip().split("\t")[3].strip()) + 3
                orf_end = int(i.strip().split("\t")[4].strip())
                orf_coordinate = chr + ":" + str(orf_start) + "-" + str(orf_end)

                # check if transcript contains annotated CDS
                if transcript in cds_coordinates.keys():
                    cds_end = int(cds_coordinates[transcript][0].split("-")[1])  # reversing coordinates for easy comparison
                    cds_start = int(cds_coordinates[transcript][-1].split("-")[0].split(":")[1])
                    if orf_start >= cds_start and orf_end <= cds_end:
                        transcript_orf_map.setdefault(transcript, []).append(protein_seq + "|" + orf_id)  # store all orfs of transcript
                else:
                    transcript_orf_map.setdefault(transcript, []).append(protein_seq + "|" + orf_id)  # store all orfs of transcript

            elif strand == "+":
                orf_start = int(i.strip().split("\t")[3].strip())
                orf_end = int(i.strip().split("\t")[4].strip()) - 3  # 3 nucleiotides substracted as ORFik counts stop codon position
                orf_coordinate = chr + ":" + str(orf_start) + "-" + str(orf_end)

                # check if transcript contains annotated CDS
                if transcript in cds_coordinates.keys():
                    cds_start = int(cds_coordinates[transcript][0].split("-")[0].split(":")[1])
                    cds_end = int(cds_coordinates[transcript][-1].split("-")[1])
                    if orf_start >= cds_start and orf_end <= cds_end:
                        transcript_orf_map.setdefault(transcript, []).append(protein_seq + "|" + orf_id)  # store all orfs of transcript
                else:
                    transcript_orf_map.setdefault(transcript, []).append(protein_seq + "|" + orf_id)  # store all orfs of transcript

            # separate annnotated and unannotated proteins
            if protein_seq in uniprot.keys():
                uniprot_accession = uniprot[protein_seq].split("|")[1]
                annotated_proteins.setdefault(protein_seq, []).append(uniprot_accession + "|" + orf_id + "|" + orf_coordinate)

            elif protein_seq in refprot.keys():
                uniprot_ids = list(uniprot.values())
                refprot_id = refprot[protein_seq]
                is_id_found = any(refprot_id in id for id in uniprot_ids)
                if is_id_found:  # avoid incorrect protein ids assinged in openprot
                    continue
                else:
                    annotated_proteins.setdefault(protein_seq, []).append(refprot_id + "|" + orf_id + "|" + orf_coordinate)
            else:
                unannotated_proteins.setdefault(protein_seq, []).append(orf_id + "|" + orf_coordinate)
                unannotated_protein_coordinates[orf_id + "|" + orf_coordinate] = protein_seq

    # find longest orf in transcript
    transcript_longest_orf_map = {}  # stores longest ORF in transcript: k: orf_id v:protein_seq
    for k in transcript_orf_map.keys():
        orfs = []
        id = []
        for orf_rec in transcript_orf_map[k]:
            orfs.append(orf_rec.split("|")[0].strip())
            id.append(orf_rec.split("|")[1].strip())
        orf_lengths = [len(orf) for orf in orfs]
        longest_orf_index = orf_lengths.index(max(orf_lengths))
        longest_orf = orfs[longest_orf_index]
        longest_orf_id = id[longest_orf_index]
        transcript_longest_orf_map[longest_orf_id] = longest_orf

    # store sequences in the temparory file which is required to perform clustering using cdhit
    orf_temp_file = os.path.join(arg_outdir, "orf_temp.txt")
    with open(orf_temp_file, 'w') as f: # Temparary file to store ORF before clustering
        entries = ''
        for protein_seq, protein_info in unannotated_proteins.items():
            entries += f">{protein_info[0]}\n" + f"{protein_seq}\n"
        for protein_seq, protein_info in annotated_proteins.items():
            entries += f">{protein_info[0]}\n" + f"{protein_seq}\n"
        _ = f.write(entries)

    # Perform cdhit clustering on unannotated proteins to consider longest representative protein sequence
    orf_clus_file = os.path.join(arg_outdir, "cdhit_out")
    SeqClust(orf_temp_file, orf_clus_file)

    orfbiotype_transcript_map = {}  # k:orf, v:transcript|biotype
    # annotations of known proteins
    orf_metadata_annotated_proteins = {}  # key:protein_seq, val:metadata
    metadata_records = []

    for k in annotated_proteins.keys():
        protein_seq = k
        transcripts = []
        gene_accession = []
        gene_symbol = []
        orf_coordinates = []
        for orf_id in annotated_proteins[k]:

            transcript = orf_id.split("|")[1].split("_")[0]
            if transcript in transcript_biotypes.keys():
                transcripts.append(transcript)
                transcript_biotype = transcript_biotypes[transcript]
                transcript_coordinates = transcript_genome_coordinates[transcript]
                gene_id = transcript_gene_id_map[transcript]
                gene_accession.append(gene_id)
                gene_name = transcript_gene_name_map[transcript]
                if re.match(r".*\S.*", gene_name.strip()):
                    gene_symbol.append(gene_name.strip())
                else:
                    gene_name = "-"  # if there is no gene name
                strand = transcript_strand[transcript]
                orf_coordinate = orf_id.split("|")[2].strip()
                orf_coordinates.append(orf_coordinate)  # orf genome coordinates
                longest_orf = "N"
                if orf_id.split("|")[1] in transcript_longest_orf_map.keys():
                    longest_orf = "Y"

                # uniprot proteins
                if protein_seq in uniprot.keys():

                    if uniprot[protein_seq].split("|")[0] == "sp":
                        protein_status = "reviewed(Swiss-Prot)"
                    elif uniprot[protein_seq].split("|")[0] == "tr":
                        protein_status = "unreviewed(TrEMBL)"

                    protein_accession = uniprot[protein_seq].split("|")[1]
                    protein_description = uniprot[protein_seq].split("|")[2]

                    orfbiotype_transcript_map.setdefault(protein_accession, []).append(transcript + "|CDS")
                    # annotate ORFs

                    get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "annotated", "CDS", "-", longest_orf, protein_status, orf_metadata_annotated_proteins, reading_frame_map[transcript], var_ORF_anno)

                # Refseq/Ensembl proteins
                elif protein_seq in refprot.keys():
                    protein_status = "-"
                    protein_accession = refprot[protein_seq]
                    protein_description = "-"
                    orfbiotype_transcript_map.setdefault(protein_accession, []).append(transcript + "|CDS")

                    # annotate ORFs
                    get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "annotated", "CDS", "-", longest_orf, protein_status, orf_metadata_annotated_proteins, reading_frame_map[transcript], var_ORF_anno)
    # annotated_proteins close

    # write output to metadata and fasta file
    canonical_annotated_orf_map = {}  # k:protein_seq, v:fasta header annotations
    all_annotated_orf_map = {}  # k:protein_seq, v:fasta header annotations

    for protein_seq, annotations in orf_metadata_annotated_proteins.items():
        for protein_anno in annotations:
            # from metadata

            protein_accession = protein_anno.split("\t")[0]
            gene_id = protein_anno.split("\t")[1]
            gene_name = protein_anno.split("\t")[2]
            protein_description = protein_anno.split("\t")[3]
            transcript = protein_anno.split("\t")[4]
            orf_coordinate = protein_anno.split("\t")[8]
            longest_orf = protein_anno.split("\t")[14]

            if arg_canonical_or_all == "canonical":
                if longest_orf == "Y":
                    metadata_records.append(protein_anno)
                    fa_header = protein_accession + "|" + orf_coordinate + "|" + gene_id + "|" + gene_name + "|" + transcript + "|" + protein_description
                    canonical_annotated_orf_map.setdefault(protein_seq, []).append(fa_header)
            elif arg_canonical_or_all == "all":
                metadata_records.append(protein_anno)
                fa_header = protein_accession + "|" + orf_coordinate + "|" + gene_id + "|" + gene_name + "|" + transcript + "|" + protein_description
                all_annotated_orf_map.setdefault(protein_seq, []).append(fa_header)

    if arg_canonical_or_all == "canonical":
        # generate fasta database
        formatFastaheader(canonical_annotated_orf_map, fw_proteomedb)
    elif arg_canonical_or_all == "all":
        # generate fasta database
        formatFastaheader(all_annotated_orf_map, fw_proteomedb)

    # Novel proteins
    AN = Annotations()
    orf_annotation_map = {}  # stores ORF annotations k:temparory ORF_id, annotations
    counter = 1  # counter for temparary ORF ids
    protein_description = "-"

    cdhit_clustering_output_file = os.path.join(arg_outdir, "cdhit_out.clstr")  # cluster file
    with open(cdhit_clustering_output_file, 'r') as cdhit_clustering_output:
        for i in cdhit_clustering_output:
            if "*" in i.strip():
                if ">Bambu" in i.strip() or ">ENS" in i.strip():
                    pattern = r">ENS.*P.*"
                    if re.search(pattern, i.strip()):  # remove representative sequences of known proteins
                        continue
                    else:
                        longest_seq_orf_id = i.strip().split(">")[1].split("...")[0]  # orf_id of Longest representative ORF in cluster
                        protein_status = "-"
                        protein_accession = "ORF_" + str(counter)
                        protein_seq = unannotated_protein_coordinates[longest_seq_orf_id]  # Protein sequence
                        # check if protein annotated in openprot
                        if protein_seq in openprot.keys():
                            openprot_id = openprot[protein_seq]
                        else:
                            openprot_id = "-"

                        for orf_id in unannotated_proteins[protein_seq]:  # accessing all ORFs coordinates for a given protein sequence
                            transcript = orf_id.split("|")[0].split("_")[0]
                            if transcript in transcript_biotypes.keys():
                                transcript_biotype = transcript_biotypes[transcript]
                                transcript_coordinates = transcript_genome_coordinates[transcript]
                                gene_id = transcript_gene_id_map[transcript]

                                if re.match(r".*\S.*", transcript_gene_name_map[transcript].strip()):
                                    gene_name = transcript_gene_name_map[transcript].strip()
                                else:
                                    gene_name = "-"

                                longest_orf = "N"
                                if orf_id.split("|")[0] in transcript_longest_orf_map.keys():
                                    longest_orf = "Y"

                                strand = transcript_strand[transcript]
                                orf_coordinate = orf_id.split("|")[1].strip()

                                if transcript_biotype == "protein_coding" and transcript in utr_coordinates.keys():

                                    utr_orf = AN.UTRAnnotations(transcript, utr_coordinates[transcript], orf_coordinate, strand, cds_coordinates[transcript])
                                    if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf:  # ORF overlap with UTR region
                                        if len(protein_seq) < arg_orf_length:   # for short uORFs, remove CDS annotation
                                            utr_orf = utr_orf.split(":")[0]
                                            # annotate ORF
                                            get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", utr_orf, openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)
                                        else:
                                            # annotate ORF
                                            get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", utr_orf, openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)

                                elif transcript_biotype != "protein_coding":

                                    if transcript in utr_coordinates.keys():  # exceptions: some transcripts have UTR but they are not protein coding

                                        utr_orf = AN.UTRAnnotations(transcript, utr_coordinates[transcript], orf_coordinate, strand, cds_coordinates[transcript])
                                        if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf:  # ORF overlap with UTR region
                                            if len(protein_seq) < arg_orf_length:
                                                utr_orf = utr_orf.split(":")[0]
                                                # annotate ORF
                                                get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", utr_orf, openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)
                                            else:
                                                # annotate ORF
                                                get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", utr_orf, openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)

                                        else:  # if ORF is not in UTR regions
                                            gene_overlap = AN.isIntergenic(orf_coordinate, protein_coding_gene_coordinates)

                                            if gene_overlap == "gene_overlap":
                                                # annotate ORF
                                                get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", "gene_overlap", openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)
                                            else:
                                                # annotate ORF
                                                get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", "intergenic", openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)

                                    # noncoding RNA transcripts. transcript doesn't have UTRs
                                    else:
                                        gene_overlap = AN.isIntergenic(orf_coordinate, protein_coding_gene_coordinates)
                                        if gene_overlap == "gene_overlap":
                                            # annotate ORF
                                            get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", "gene_overlap", openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)
                                        else:
                                            # annotate ORF
                                            get_protein_annotation(transcript, gene_id, gene_name, protein_description, var_transcript_ORF_map, protein_seq, protein_accession, strand, transcript_biotype, transcript_coordinates, orf_coordinate, "unannotated", "intergenic", openprot_id, longest_orf, protein_status, orf_annotation_map, reading_frame_map[transcript], var_ORF_anno)

                        # for loop closed
                        counter = counter + 1

    # write output of novel ORFs
    canonical_novel_orf_ids = []
    all_novel_orf_ids = []

    # store temparary ids of canonical or all novel ORFs according to options provided
    for protein_seq, annotations in orf_annotation_map.items():
        for protein_annotation in annotations:
            longest_orf = protein_annotation.split("\t")[14]
            accession = protein_annotation.split("\t")[0].replace("_var", "")
            if arg_canonical_or_all == "canonical":
                if longest_orf == "Y":
                    if accession not in canonical_novel_orf_ids:
                        canonical_novel_orf_ids.append(accession)
            if arg_canonical_or_all == "all":
                if accession not in all_novel_orf_ids:
                    all_novel_orf_ids.append(accession)
    # loop closed

    # replace the temporary ORFs ids based on their index
    canonical_novel_orf_map = {}  # k:protein_seq, v:fasta header annotations
    all_novel_orf_map = {}  # k:protein_seq, v:fasta header annotations

    for protein_seq, annotations in orf_annotation_map.items():
        for protein_annotation in annotations:
            accession = protein_annotation.split("\t")[0]
            gene_id = protein_annotation.split("\t")[1]
            gene_name = protein_annotation.split("\t")[2]
            protein_description = protein_annotation.split("\t")[3]
            transcript = protein_annotation.split("\t")[4]
            orf_coordinate = protein_annotation.split("\t")[8]
            longest_orf = protein_annotation.split("\t")[14]
            if arg_canonical_or_all == "canonical":
                if longest_orf == "Y":
                    if "_var" in accession:
                        temp_acc = accession.replace("_var", "")
                        new_accesion = "ORF_" + str(canonical_novel_orf_ids.index(temp_acc) + 1) + "_var"  # +1 as list index starts with zero
                    else:
                        new_accesion = "ORF_" + str(canonical_novel_orf_ids.index(accession) + 1)  # +1 as list index starts with zero

                    fa_header = new_accesion + "|" + orf_coordinate + "|" + gene_id + "|" + gene_name + "|" + transcript

                    canonical_novel_orf_map.setdefault(protein_seq, []).append(fa_header)
                    revised_annotations = protein_annotation.replace(accession, new_accesion)
                    metadata_records.append(revised_annotations)

            if arg_canonical_or_all == "all":
                if "_var" in accession:
                    temp_acc = accession.replace("_var", "")
                    new_accesion = "ORF_" + str(all_novel_orf_ids.index(temp_acc) + 1) + "_var"  # +1 as list index starts with zero
                else:
                    new_accesion = "ORF_" + str(all_novel_orf_ids.index(accession) + 1)

                fa_header = new_accesion + "|" + orf_coordinate + "|" + gene_id + "|" + gene_name + "|" + transcript
                all_novel_orf_map.setdefault(protein_seq, []).append(fa_header)
                revised_annotations = protein_annotation.replace(accession, new_accesion)
                metadata_records.append(revised_annotations)
    # loop closed

    metadata_records_uniq = list(dict.fromkeys(metadata_records))

    # write metadata
    proteome_database_metadata_filename = os.path.join(arg_outdir, "proteome_database_metadata.txt")
    with open(proteome_database_metadata_filename, 'w') as f:
        f.write('\t'.join(["accession", "gene", "gene_symbol", "protein_description", "transcript", "strand", "transcript_biotype", "transcript_coordinates", "orf_genomic_coordinates", "reading_frame", "orf_type", "localisation", "openprot_id", "protein_sequence", "longest_orf_in_transcript", "uniprot_status", "amino_acid_change", "molecular_weight(kDA)", "isoelectric_point", "hydrophobicity", "aliphatic_index"]) + '\n')
        f.writelines(metadata_records_uniq)

    # write fasta output
    if arg_canonical_or_all == "canonical":
        # generate fasta database
        formatFastaheader(canonical_novel_orf_map, fw_proteomedb)
    elif arg_canonical_or_all == "all":
        # generate fasta database
        formatFastaheader(all_novel_orf_map, fw_proteomedb)

    # remove intermediate files
    remove_file_if_exists(orf_temp_file)
    remove_file_if_exists(orf_clus_file)
    remove_file_if_exists(cdhit_clustering_output_file)

if __name__ == "__main__":
    main()