from pycdhit import cd_hit, read_clstr
import peptides as pep
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

def SeqClust(fasta, out):
    cd_hit(i = fasta, o = out, c = 1.0, d = 0, sc = 1, n = 3, G = 0, aS = 1, aL = 0, A = 0)

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

class Annotations():
    def UTRAnnotations(self,tr,utrmap,orf,st,trcds): #transcript,utrcord, orfcoord,strand,transcript coordinates
        if ":" in orf and "-" in orf:
          orfBT=""
          ochr=orf.split(':')[0]
          utr_anno={} #annotates 5' and 3' UTR coordinates: k: UTR coordinates v:3UTR/5UTR
          utr_orf_pos=[]
          #segreagate 3UTR and 5UTR coordinates and determine postion of ORFs in UTR exons
          if st=="+":
            ostart=int(orf.split(":")[1].split('-')[0])
            oend=int(orf.split(":")[1].split('-')[1])
            first_cds=int(trcds[0].split(':')[1].split("-")[0]) #start position of first CDs of transcript
            last_cds=int(trcds[-1].split(":")[1].split("-")[1]) # end position of last CDS of transcript
            for utr in utrmap:  #Identifying 5' UTR and 3'UTR exons among UTR cooordinates of transcript
              uchr=utr.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utr.split(':')[1].split('-')[0])
                u_end=int(utr.split(':')[1].split('-')[1])
                #marking 5UTR and 3UTR
                if u_start < first_cds and u_end < first_cds:#if UTR start and end is lesser than first cds start then it is 5UTR
                  utr_anno[utr]='5UTR'
                elif u_start > last_cds and u_end>last_cds: #if UTR start and end is greater than first cds then it is 5UTR
                  utr_anno[utr]='3UTR'

            for utrs in utrmap:
              uchr=utrs.split(':')[0].strip()
              u_start=int(utrs.split(':')[1].split('-')[0])
              u_end=int(utrs.split(':')[1].split('-')[1])
              if uchr==ochr:
                if ostart >=u_start and ostart < u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("*"+utr_anno[utrs]) # * indicates if ORF starts in 5' or 3' UTR
                elif ostart >u_start and ostart <= u_end:
                  utr_orf_pos.append("*"+utr_anno[utrs])

                if oend >= u_start and oend< u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR
                elif oend > u_start and oend<= u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR

          if st=="-":
            ostart=int(orf.split(":")[1].split('-')[1])
            oend=int(orf.split(":")[1].split('-')[0])
            first_cds=int(trcds[0].split(':')[1].split("-")[1])
            last_cds=int(trcds[-1].split(":")[1].split("-")[0])
            for utr in utrmap:
              uchr=utr.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utr.split(':')[1].split('-')[1])
                u_end=int(utr.split(':')[1].split('-')[0])
                if u_start > first_cds and u_end > first_cds :
                  utr_anno[utr]='5UTR'
                elif u_start < last_cds and u_end<last_cds:
                  utr_anno[utr]='3UTR'

            for utrs in utrmap:
              uchr=utrs.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utrs.split(':')[1].split('-')[1])
                u_end=int(utrs.split(':')[1].split('-')[0])
                #start
                if ostart <= u_start and ostart > u_end: #if ORF start is lesser than UTR start and greater than than UTR end
                  utr_orf_pos.append("*"+utr_anno[utrs]) # * indicates if ORF starts in 5' or 3' UTR
                elif ostart < u_start and ostart >= u_end:
                  utr_orf_pos.append("*"+utr_anno[utrs])
                #end
                if oend<= u_start and oend > u_end: #if ORF end is lesser than UTR start and greater than than UTR end
                    utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR
                elif oend< u_start and oend >= u_end:
                  utr_orf_pos.append("**"+utr_anno[utrs])

          # Determine ORF type
          if len(utr_orf_pos)==1:
            if utr_orf_pos[0]=='*5UTR':
              orfBT="5UTR:CDS"
            elif utr_orf_pos[0]=='*3UTR':
              orfBT="3UTR"
            elif utr_orf_pos[0]=='**3UTR':
              orfBT="CDS:3UTR" #Alt-CDS:3UTR
            elif utr_orf_pos[0]=='**5UTR': #handeling ambiguity due to ORFs locations predicted by ORFik are off by certain nucleiotide differences
              orfBT="5UTR"

          elif len(utr_orf_pos)==2:
            if utr_orf_pos[0]=='*3UTR' and utr_orf_pos[1]=='**3UTR':
              orfBT="3UTR"
            elif utr_orf_pos[0]=='*5UTR' and utr_orf_pos[1]=='**5UTR':
              orfBT="5UTR"
            elif utr_orf_pos[0]=='*5UTR' and utr_orf_pos[1]=='**3UTR':
              orfBT="5UTR:3UTR"
            elif utr_orf_pos[0]=='**3UTR' and utr_orf_pos[1]=='*5UTR':
              orfBT="5UTR:3UTR"
          return orfBT

    def isIntergenic(self,orf_cord,gene_cord_map):
      #ORF
      if ":" in orf_cord and "-" in orf_cord:
        ochr=orf_cord.split(":")[0]
        ostart=int(orf_cord.split(":")[1].split("-")[0])
        oend=int(orf_cord.split(":")[1].split("-")[1])
        orf_overlap=""
        for gnpost in gene_cord_map.values():
          #gene
          genchr=gnpost.split(":")[0]
          gnstart=int(gnpost.split(":")[1].split('-')[0])
          gnend=int(gnpost.split(":")[1].split('-')[1])
          if ochr==genchr:
            if ostart >= gnstart and oend <= gnend:
              orf_overlap="gene_overlap"
              break
        return orf_overlap

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

def calculate_similarity(seq1, seq2):
    # Use a substitution matrix like BLOSUM62
    matrix = matlist.blosum62

    # Set gap open and gap extend penalties
    gap_open = -10  # Penalty for opening a gap
    gap_extend = -0.5  # Penalty for extending a gap

    # Perform global alignment with the matrix and gap penalties
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    # Take the best alignment (the first one)
    best_alignment = alignments[0]

    # Extract aligned sequences
    aligned_seq1, aligned_seq2, score, start, end = best_alignment

    # Calculate the percentage similarity
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    similarity_percentage = matches / len(aligned_seq1) * 100

    # Find changed amino acids
    changed_positions = [(i + 1, a, b) for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)) if a != b and a != "-" and b != "-"]  # first mutant and second wild type

    # convert into standard format to denote mutation
    aa_changes = ",".join([f"{b}{pos}{a}" for pos, a, b in changed_positions])

    return similarity_percentage, aa_changes