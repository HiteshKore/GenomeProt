from pycdhit import cd_hit
import peptides as pep
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

def SeqClust(fasta, out, num_threads, memory_limit):
    cd_hit(i = fasta, o = out, c = 1.0, d = 0, sc = 1, n = 3, G = 0, aS = 1, aL = 0, A = 0, T = num_threads, M = memory_limit)

def getDbAnnotations(db_filename, openprot_map, ref_prot_map, uniprot_map):  # Reference protein sequence database
    with open(db_filename, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if '#' in line:
                continue

            columns = [col.strip() for col in line.split('\t')]
            [protein_id, seq, uniprot_header] = columns[:3]
            protein_id = protein_id.split('|')[0].replace('>', '')

            openprot_map[seq] = protein_id
            if ">IP_" not in line and ">II_" not in line:
                ref_prot_map[seq] = protein_id
            if uniprot_header != '-':
                uniprot_accession = uniprot_header.split('|')[0].replace('>', '') + '|' + uniprot_header.split('|')[1]  # trEMBL or Reviewed| Accession
                uniprot_gene_description = columns[3]
                uniprot_map[seq] = f"{uniprot_accession}|{uniprot_gene_description}"

def UTRAnnotations(utrmap, orf, st, trcds): # utrcoord, orfcoord, strand, transcript coordinates
    if not (':' in orf and '-' in orf and (st == '+' or st == '-')):
        return ''

    is_forward_strand = (st == '+')
    [ochr, orf_info] = [col.strip() for col in orf.split(':')[:2]]
    [orf_start, orf_end] = [int(col) for col in orf_info.split('-')[:2]]
    if not is_forward_strand:
        orf_start, orf_end = orf_end, orf_start

    first_cds = int(trcds[0].split(':')[1].split('-')[0 if is_forward_strand else 1])  # start position of first CDS of transcript
    last_cds = int(trcds[-1].split(':')[1].split('-')[1 if is_forward_strand else 0])  # end position of last CDS of transcript

    # segreagate 3UTR and 5UTR coordinates and determine postion of ORFs in UTR exons
    utr_orf_pos = []
    for utr in utrmap:  # Identifying 5' UTR and 3'UTR exons among UTR cooordinates of transcript
        uchr = utr.split(':')[0].strip()
        if uchr != ochr:
            continue

        [utr_start, utr_end] = [int(col) for col in utr.split(':')[1].split('-')[:2]]
        if not is_forward_strand:
            utr_start, utr_end = utr_end, utr_start

        utr_annotation = ''
        if (    is_forward_strand and utr_start < first_cds and utr_end < first_cds) or \
           (not is_forward_strand and utr_start > first_cds and utr_end > first_cds):
            utr_annotation = "5UTR"
        elif (    is_forward_strand and utr_start > last_cds and utr_end > last_cds) or \
             (not is_forward_strand and utr_start < last_cds and utr_end < last_cds):
            utr_annotation = "3UTR"
        if utr_annotation == '':
            continue

        if (    is_forward_strand and ((utr_start <= orf_start < utr_end) or (utr_start < orf_start <= utr_end))) or \
           (not is_forward_strand and ((orf_start <= utr_start and orf_start > utr_end) or (orf_start < utr_start and orf_start >= utr_end))):
            utr_orf_pos.append('*' + utr_annotation)
        if (    is_forward_strand and ((utr_start <= orf_end < utr_end) or (utr_start < orf_end <= utr_end))) or \
           (not is_forward_strand and ((orf_end <= utr_start and orf_end > utr_end) or (orf_end < utr_start and orf_end >= utr_end))):
            utr_orf_pos.append("**" + utr_annotation)

    # Determine the ORF type
    orfBT = ""
    if utr_orf_pos == ["*5UTR"]:
        orfBT = "5UTR:CDS"
    elif utr_orf_pos == ["**3UTR"]:
        orfBT = "CDS:3UTR"  # Alt-CDS:3UTR
    elif utr_orf_pos == ["**5UTR"] or utr_orf_pos == ["*5UTR", "**5UTR"]:   # handling ambiguity due to ORF locations predicted by ORFik being off by certain nucleotide differences
        orfBT = "5UTR"
    elif utr_orf_pos == ["*3UTR"] or utr_orf_pos == ["*3UTR", "**3UTR"]:
        orfBT = "3UTR"
    elif utr_orf_pos == ["*5UTR", "**3UTR"] or utr_orf_pos == ["**3UTR", "*5UTR"]:
        orfBT = "5UTR:3UTR"
    return orfBT

def isIntergenic(orf_coord, gene_coord_map):
    [orf_chrom, orf_range] = orf_coord.split(':')[:2]
    [orf_start, orf_end] = [int(num) for num in orf_range.split('-')]
    for (gene_start, gene_end) in gene_coord_map.get(orf_chrom, []):
        if gene_start <= orf_start and orf_end <= gene_end:
            return True
    return False

# function to calculate protein sequence physico-chemical properties
def calculate_sequence_properties(protein_seq):
    peptide_obj = pep.Peptide(protein_seq)
    Mol_wt = calculateMolWt(peptide_obj)
    IsoPt = calculateIsoElectricPoint(peptide_obj)
    HydroP_ind = calculateHydrophobicity(peptide_obj)
    Aliphatic_ind = calculateAliphatic_index(peptide_obj)
    seq_properties = '\t'.join(map(str, [Mol_wt, IsoPt, HydroP_ind, Aliphatic_ind]))
    return seq_properties

def calculateMolWt(peptide_obj):
    return round(float(peptide_obj.molecular_weight()) / 1000, 2)   # in KDa

def calculateIsoElectricPoint(peptide_obj):
    return round(peptide_obj.isoelectric_point(), 2)

def calculateHydrophobicity(peptide_obj):
    return round(peptide_obj.hydrophobicity(), 2)

def calculateAliphatic_index(peptide_obj):
    return round(peptide_obj.aliphatic_index(), 2)

matrix = MatrixInfo.blosum62
def calculate_similarity(seq1, seq2):
    # Perform global alignment with the BLOSUM62 matrix, -10 gap open penalty and -0.5 gap extend penalty
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)

    # Take the best alignment (the first one)
    best_alignment = alignments[0]

    # Extract aligned sequences
    aligned_seq1, aligned_seq2 = best_alignment[:2]

    # Calculate the percentage similarity
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    similarity_percentage = matches / len(aligned_seq1) * 100

    # Find changed amino acids
    changed_positions = [(i + 1, a, b) for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)) if a != b and a != "-" and b != "-"]  # first mutant and second wild type

    # convert into standard format to denote mutation
    aa_changes = ",".join([f"{b}{pos}{a}" for pos, a, b in changed_positions])

    return similarity_percentage, aa_changes