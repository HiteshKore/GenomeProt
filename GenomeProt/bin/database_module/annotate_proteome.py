#This script removes the redundant ORFs and considers longest ORF if ORF is part of longer ORF.
#It annotates them based on their location on the genome/transcript
#cd-hit (First install cd-hit commond line tool in linux (command to install:'conda install bioconda/label/cf201901::cd-hit'). Next, install py-cdhit library using  'pip install py-cdhit' command)
#Usage:python3 annotate_proteome.py gencode.vM33.chr_patch_hapl_scaff.annotation_chrX.gtf openprot_uniprotDb_mm.txt  ORFome_aa.txt proteome_database_transcripts.gtf <outdir> <canonical/all> <orf_length> <M> <organism:HUMAN,CAEEL,MOUSE,RAT,DROME,DANRE>
#######################################################################

import sys
import os
import re
import peptides as pep
from parse_reference_gtf  import *
from annotate_proteome_functions import *
import subprocess


def main():
  args=sys.argv
  if len(args) != 10:
        print("Usage: python3 annotate_proteome.py <reference_gtf> <custom_openprot+uniprot_db> <ORFome_aa.txt> <ORFome_transcripts.gtf> <outdir> <canonical/all> <orf_length> <variant_protein_db/None> <organism>")
        sys.exit(1)
  else:
    
    #custom openprot+ uniprot annotation database
    RF=refDb()
    uniprot={} #comprise UniProt annotations: k:seq v: trEMBL/reviewed|protein_id| gene description
    openprot={} #comprise OpenProt annotations: k:seq v:protein_id
    refprot={} #comprise reference proteins annotated in uniprot/refseq/ensembl. k:seq v:protein id
    RF.getDbAnnotations(args[2],openprot,refprot,uniprot) #OpenProt annotations
    
    #GTF annotations
    #Reference GTF ENSEMBL/GENCODE
    refence_gtf=open(args[1])
    
    #Gene
    gene_biotype={} #k:gene_id i.e ENSG,k:gene biotype
    gene_coordinates={} #k:gene_id,k:gene coordinates
    protein_coding_gene_coordinates={} #k:gene_id,k:gene coordinates
    
    #transcripts
    transcript_biotypes={} #k:transcript_id i.e ENST,k: transcript biotype
    transcript_genome_coordinates={} #k:transcript_id i.e ENST,k: genome coordinates
    transcript_gene_id_map={} #k:transcript_id i.e ENST,k: gene id i.e. ENSG
    transcript_gene_name_map={} #k:transcript_id i.e ENST,k: gene name
    transcript_strand={} #k:transcript_id i.e ENST,k: strand
    #exon
    exon_lengths={} #transcript exon lengths k:transcript_id i.e ENST,k: lenghts of exons
    exon_coordinates={} #genome coordinates of transcript exons k:transcript_id i.e ENST,k: genomic coordinates
    
    #utr
    utr_coordinates={} #genome coordinates of transcript utrs k:transcript_id i.e ENST,k: genomic coordinates
    #cds coordinates
    cds_coordinates={} #genome coordinates of transcript cds k:transcript_id i.e ENST,k: genomic coordinates
    
    for i in refence_gtf:
      if i.startswith("#"):
        continue
      else:
        GP=GTFParser(i.strip())
        if GP.feature=="gene":
          gene_biotype[GP.ensgene]=GP.genetype
          gene_coordinates[GP.ensgene]=GP.getCoordinates()
          if GP.genetype=="protein_coding": #Feching coordinates for protein coding genes
            protein_coding_gene_coordinates[GP.ensgene]=GP.getCoordinates()
        if GP.feature=="transcript":
          transcript_biotypes[GP.transcriptid]=GP.transcript_type
          transcript_genome_coordinates[GP.transcriptid]=GP.getCoordinates()
          transcript_strand[GP.transcriptid] = GP.strand
          transcript_gene_id_map[GP.transcriptid]=GP.ensgene
          transcript_gene_name_map[GP.transcriptid]=GP.genename
        
        if GP.featureExists('exon'):
          exon_lengths.setdefault(GP.transcriptid, []).append(GP.end - GP.start+1) #stores exon length
          exon_coordinates.setdefault(GP.transcriptid,[]).append(str(GP.start)+"-"+str(GP.end)) #stores exon gen
        if GP.featureExists('UTR'):
          utr_coordinates.setdefault(GP.transcriptid,[]).append(GP.getCoordinates())
        if GP.featureExists('CDS'):
          cds_coordinates.setdefault(GP.transcriptid,[]).append(GP.getCoordinates())
    refence_gtf.close()
    
    #ORFome_transcripts gtf
    
    orfome_transcript_gtf=open(args[4])
    for i in orfome_transcript_gtf:
      if i.startswith("#"):
        continue
      else:
        if 'BambuGene' in i.strip() or 'BambuTx' in i.strip():#denovoGene denovoTx
          GP =GTFParser(i.strip())
          if GP.featureExists('exon'):
            exon_lengths.setdefault(GP.transcriptid, []).append(GP.end - GP.start+1) #stores exon length
            exon_coordinates.setdefault(GP.transcriptid,[]).append(str(GP.start)+"-"+str(GP.end)) #stores exon gen
          if GP.featureExists('transcript'):
            transcript_genome_coordinates[GP.transcriptid] = GP.getCoordinates()
            transcript_strand[GP.transcriptid] = GP.strand
            transcript_gene_id_map[GP.transcriptid]=GP.ensgene
            transcript_gene_name_map[GP.transcriptid]=GP.genename
            if 'BambuTx' in GP.transcriptid:
              transcript_biotypes[GP.transcriptid]="novel"
            if 'BambuGene' in GP.genename:
              gene_biotype[GP.genename]="novel"
            
    orfome_transcript_gtf.close()
    
    
    var_transcript_ORF_map={}
    var_ORF_anno={}
    if args[8] !="None" and os.path.isfile(args[8]): #variant file is optional
      mutant_protein_fh=open(args[8]) 
      for i in mutant_protein_fh:
        if i.strip().startswith("transcript"):
          continue
        else:
          transcript=i.strip().split("\t")[0]
          if transcript in transcript_gene_id_map.keys(): #the gencode version of reference GTF must be same
            gene_id=transcript_gene_id_map[transcript]
            var_sq=i.strip().split("\t")[1]
            orf_coordinate=i.strip().split("\t")[2]
            mutation_type=i.strip().split("\t")[3]
            var_transcript_ORF_map.setdefault(transcript,[]).append(var_sq)
            var_ORF_anno[var_sq]=orf_coordinate+"|"+mutation_type
      
      mutant_protein_fh.close()
    
    #function to calculate protein sequence physico-chemical properties
    def calculate_sequence_properties(protein_seq):
      SP=SequenceProperties()
      Mol_wt=SP.calculateMolWt(protein_seq)
      IsoPt=SP.calculateIsoElectricPoint(protein_seq)
      HydroP_ind=SP.calculateHydrophobicity(protein_seq)
      Aliphatic_ind=SP.calculateAliphatic_index(protein_seq)
      seq_properties=str(Mol_wt)+"\t"+str(IsoPt)+"\t"+str(HydroP_ind)+"\t"+str(Aliphatic_ind)
      return seq_properties
    
    #function to identify amino acid residue changes in wildtype and variant protein sequence
    def variant_protein_annotations(transcript,var_orfs,wt_protein):
      res=""
      for var_protein in var_orfs:
        SS=SequenceSimilarity()
        similarity, aa_change = SS.calculate_similarity(var_protein, wt_protein) #parise wise alignment to calculate sequence similarity
        similarity=round(similarity,2)
        if similarity >96 and similarity <100:
          sq_properties=calculate_sequence_properties(var_protein)
          coordinates=var_ORF_anno[var_protein].split("|")[0]
          var_type=var_ORF_anno[var_protein].split("|")[1]
          res =var_protein+"\t"+aa_change+"\t"+coordinates+"\t"+var_type+"\t"+sq_properties
      return res
    
    #function to annotate protein sequence based on genomic coordinates, biotypes, physioco-chemical properties
    def get_protein_annotation(transcript,gene_id,gene_name,protein_des,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,orf_type,localisation,openprot_annotations,longest_orf,protein_status,orf_metadata_map):
      if strand=="-":
        orf_start=int(orf_coordinate.strip().split(':')[1].split("-")[1])
        rf=(orf_start - 1) % 3
      else:
        orf_start=int(orf_coordinate.strip().split(':')[1].split("-")[0])
        rf=(orf_start - 1) % 3
      
      if transcript in var_transcript_ORF_map.keys():
        var_prot_attr=variant_protein_annotations(transcript,var_transcript_ORF_map[transcript],protein_seq) #identify variants and calculate physico-chemical properties
        if var_prot_attr!="":
          var_seq=var_prot_attr.split("\t")[0]
          variants=var_prot_attr.split("\t")[1]
          var_seq_coordinates=var_prot_attr.split("\t")[2]
          chr=var_seq_coordinates.strip().split(':')[0]
          if strand=="-":
            orf_start=int(var_seq_coordinates.strip().split(':')[1].split("-")[0])+3
            orf_end=var_seq_coordinates.strip().split(':')[1].split("-")[1]
            rf=(int(orf_end) - 1) % 3
            var_seq_coordinates=chr+":"+str(orf_start)+"-"+orf_end
            print(protein_accession+"_var",var_seq_coordinates,strand,orf_end)
          elif strand=="+":
            orf_start=var_seq_coordinates.strip().split(':')[1].split("-")[0]
            orf_end=int(var_seq_coordinates.strip().split(':')[1].split("-")[1])-3 # 3 nucleiotides substracted as ORFik counts stop codon position
            rf=(int(orf_start) - 1) % 3
            var_seq_coordinates=chr+":"+orf_start+"-"+str(orf_end)
            print("*****",protein_accession+"_var",var_seq_coordinates,strand)

          var_type=var_prot_attr.split("\t")[3] #variants
          var_seq_prop=var_prot_attr.split("\t")[4] #sequence properties
          if var_type=="HM": #homozygous variants
            orf_annotation=protein_accession+"_var"+"\t"+gene_id+"\t"+gene_name+"\t"+protein_des+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+var_seq_coordinates+"\t"+str(rf)+"\t"+orf_type+"\t"+localisation+"\t"+openprot_annotations+"\t"+var_seq+"\t"+longest_orf+"\t"+protein_status+"\t"+variants+"\t"+calculate_sequence_properties(var_seq)+"\n"
            orf_metadata_map.setdefault(var_seq,[]).append(orf_annotation)
          elif var_type=="HT": #heterozygous variants
            orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+protein_des+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\t"+str(rf)+"\t"+orf_type+"\t"+localisation+"\t"+openprot_annotations+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+calculate_sequence_properties(protein_seq)+"\n"
            orf_metadata_map.setdefault(protein_seq,[]).append(orf_annotation) #wildtype sequence
            orf_annotation=protein_accession+"_var"+"\t"+gene_id+"\t"+gene_name+"\t"+protein_des+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+var_seq_coordinates+"\t"+str(rf)+"\t"+orf_type+"\t"+localisation+"\t"+openprot_annotations+"\t"+var_seq+"\t"+longest_orf+"\t"+protein_status+"\t"+variants+"\t"+calculate_sequence_properties(var_seq)+"\n"
            orf_metadata_map.setdefault(var_seq,[]).append(orf_annotation) #variant sequence
           
      else: #if transcript not in var_transcript_ORF_map
        orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+protein_des+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\t"+str(rf)+"\t"+orf_type+"\t"+localisation+"\t"+openprot_annotations+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+calculate_sequence_properties(protein_seq)+"\n"
        orf_metadata_map.setdefault(protein_seq,[]).append(orf_annotation)
        
    #function to generate fasta sequence for proteome database
    def formatFastaheader(orf_info_map,outfile): #takes annotated ORFs map and output fasta file handler
      for prot_sq,annotations in orf_info_map.items():
        accession=""
        org=args[9]
        orgmap = {'HUMAN': 'Homo sapiens',
              'CAEEL': 'Caenorhabditis elegans',
              'MOUSE': 'Mus musculus',
              'RAT': 'Rattus norvegicus',
              'DROME': 'Drosophila melanogaster',
              'DANRE': 'Danio rerio'}
      

        orf_coord=[] #ORF coordinates
        GA=[] #gene accession
        GN=[] #gene name
        TA=[] #transcript\
        description=[]
        #protein_accession+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript+"|"+protein_description
        for header_info in annotations:
          accession=header_info.split("|")[0] #accession
          orf_coord.append(header_info.split("|")[1]) #ORF coordinates
          GA.append(header_info.split("|")[2]) #gene accession
          GN.append(header_info.split("|")[3]) #gene name
          TA.append(header_info.split("|")[4]) #transcript
          if len(header_info.split("|"))==6:
            description.append(header_info.split("|")[5])
        orf_coord_s="CO="+",".join(list(set(orf_coord)))
        TA_s="TA="+",".join(list(set(TA)))
        GA_s="GA="+",".join(list(set(GA)))
        des_s=",".join(list(set(description)))

        if len(GN)==0: #varify
          GN.append("-")
        GN_s="GN="+",".join(list(set(GN)))
        if len(list(set(description)))>0:
          fa_seq=">kn|"+accession+"|"+accession+"_"+org+" "+des_s+" OS="+orgmap[org]+" "+orf_coord_s+" "+GA_s+" "+GN_s+" "+TA_s+"\n"+prot_sq+"\n"
        else:
          if 'ORF_' in accession:
            fa_seq=">nv|"+accession+"|"+accession+"_"+org+" -"+" OS="+orgmap[org]+" "+orf_coord_s+" "+GA_s+" "+GN_s+" "+TA_s+"\n"+prot_sq+"\n"
          else:
            fa_seq=">kn|"+accession+"|"+accession+"_"+org+" -"+" OS="+orgmap[org]+" "+orf_coord_s+" "+GA_s+" "+GN_s+" "+TA_s+"\n"+prot_sq+"\n"
            
        outfile.write(fa_seq)
        
        
    #ORFomedb
    unannotated_proteins={} #stores unannotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    unannotated_protein_coordinates={} #stores unannotated protein coordinates as key: k: transcript|orf_number|genomic_coordinates,v:protein sequence
    orfome=open(args[3]) 
    transcript_orf_map={} #stores all orfs for each transcript: k;transcript_id,v:protein_seq|orf_id
    
    annotated_proteins={} #store annotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    
    for i in orfome:
      if i.startswith("ORF_id"):
        continue
      else:
        orf_id=i.strip().split('\t')[0].strip()
        protein_seq=i.strip().split('\t')[1].strip().replace("*","")
        transcript=orf_id.split("_")[0]
        
        strand=i.strip().split('\t')[5].strip()
        chr=i.strip().split('\t')[2].strip()
        if strand =="-":
          orf_start=int(i.strip().split('\t')[3].strip())+3
          orf_end=int(i.strip().split('\t')[4].strip())
          orf_coordinate=chr+":"+str(orf_start)+"-"+str(orf_end)
          
          #check if transcript contains annotated CDS
          if transcript in cds_coordinates.keys():
            #print(transcript,cds_coordinates[transcript])
            cds_end=int(cds_coordinates[transcript][0].split("-")[1]) #reversing coordinates for easy comparison
            cds_start=int(cds_coordinates[transcript][-1].split("-")[0].split(":")[1])
            if orf_start >= cds_start and orf_end <=cds_end:
              #print("orf-coord",orf_start,orf_end)
              transcript_orf_map.setdefault(transcript,[]).append(protein_seq+"|"+orf_id) #store all orfs of transcript
          else:
            transcript_orf_map.setdefault(transcript,[]).append(protein_seq+"|"+orf_id) #store all orfs of transcript

        elif strand =="+":
          orf_start=int(i.strip().split('\t')[3].strip())
          orf_end=int(i.strip().split('\t')[4].strip())-3 # 3 nucleiotides substracted as ORFik counts stop codon position
          orf_coordinate=chr+":"+str(orf_start)+"-"+str(orf_end)
          
          #check if transcript contains annotated CDS
          if transcript in cds_coordinates.keys():
            
            cds_start=int(cds_coordinates[transcript][0].split("-")[0].split(":")[1])
            cds_end=int(cds_coordinates[transcript][-1].split("-")[1])
            #print(transcript,cds_coordinates[transcript])
            if orf_start >= cds_start and orf_end <=cds_end:
              #print(cds_start,cds_end)
              #print("orf-coord",orf_start,orf_end)
              transcript_orf_map.setdefault(transcript,[]).append(protein_seq+"|"+orf_id) #store all orfs of transcript
              
          else:
            transcript_orf_map.setdefault(transcript,[]).append(protein_seq+"|"+orf_id) #store all orfs of transcript

        #separate annnotated and unannotated proteins
        if protein_seq in uniprot.keys():
          uniprot_accession=uniprot[protein_seq].split("|")[1]
          annotated_proteins.setdefault(protein_seq,[]).append(uniprot_accession+"|"+orf_id+"|"+orf_coordinate)
          
        elif protein_seq in refprot.keys():
          uniprot_ids=list(uniprot.values())
          refprot_id=refprot[protein_seq]
          is_id_found = any( refprot_id in id for id in uniprot_ids)
          if is_id_found:  #avoid incorrect protein ids assinged in openprot
            continue
          else:
            annotated_proteins.setdefault(protein_seq,[]).append(refprot_id+"|"+orf_id+"|"+orf_coordinate)
            
        else:
          unannotated_proteins.setdefault(protein_seq,[]).append(orf_id+"|"+orf_coordinate)
          unannotated_protein_coordinates[orf_id+"|"+orf_coordinate]=protein_seq

    orfome.close()
    
    #find longest orf in transcript
    transcript_longest_orf_map={} #stores longest ORF in transcript: k: orf_id v:protein_seq
    for k in transcript_orf_map.keys():
      orfs=[]
      id=[]
      for orf_rec in transcript_orf_map[k]:
        orfs.append(orf_rec.split("|")[0].strip())
        id.append(orf_rec.split("|")[1].strip())
      orf_lengths = [len(orf) for orf in orfs]
      longest_orf_index=orf_lengths.index(max(orf_lengths))
      longest_orf=orfs[longest_orf_index]
      longest_orf_id=id[longest_orf_index]
      transcript_longest_orf_map[longest_orf_id]=longest_orf
      
    #store sequences in the temparory file which is required to perform clustering using cdhit
    #Ouputdir 
    
    outdir=args[5]+"/"
    orf_temp_file=outdir+"orf_temp.txt"
    fw_orf_temp=open(orf_temp_file,"w") #Temparary file to store ORF before clustering
    
    for k in unannotated_proteins.keys():
      fw_orf_temp.write(">"+unannotated_proteins[k][0]+"\n"+k+"\n")
    
    for k in annotated_proteins.keys():
      fw_orf_temp.write(">"+annotated_proteins[k][0]+"\n"+k+"\n")
    
    fw_orf_temp.close()
    
    #cdhit clusering on unannotated proteins
    orf_clus=outdir+"cdhit_out"
    CS=clusterSequences() #Perform clustering to consider longest representative protein sequence
    clusters=CS.SeqClust(orf_temp_file,orf_clus)
    #cluster file
    cdhit_clustering_output=open(orf_clus+".clstr")

    #annotations
    fw_proteomedb=open(outdir+"proteome_database.fasta","w")
    fw_proteomedb_metadata=open(outdir+"proteome_database_metadata.txt","w")
      
    # Create a tab-separated format string
    fw_proteomedb_metadata.write("accession\tgene\tgene_symbol\tprotein_description\ttranscript\tstrand\ttranscript_biotype\ttranscript_coordinates\torf_genomic_coordinates\treading_frame\torf_type\tlocalisation\topenprot_id\tprotein_sequence\tlongest_orf_in_transcript\tuniprot_status\tamino_acid_change\tmolecular_weight(kDA)\tisoelectric_point\thydrophobicity\taliphatic_index\n")
        
    longest_orf_anno={} #k:accession v:Y, N
    orfbiotype_transcript_map={} #k:orf, v:transcript|biotype
    #annotations of known proteins
    orf_metadata_annotated_proteins={} #key:protein_seq, val:metadata

    for k in annotated_proteins.keys():
      protein_seq=k
      transcripts=[]
      gene_accession=[]
      gene_symbol=[]
      orf_coordinates=[]
      
      for orf_id in  annotated_proteins[k]:
        transcript=orf_id.split("|")[1].split("_")[0]

        if transcript in transcript_biotypes.keys():
          transcripts.append(transcript)
          transcript_biotype=transcript_biotypes[transcript]
          transcript_coordinates=transcript_genome_coordinates[transcript]
          gene_id=transcript_gene_id_map[transcript]
          gene_accession.append(gene_id)
          gene_name=transcript_gene_name_map[transcript]
          if re.match(r'.*\S.*',gene_name.strip()):
              gene_symbol.append(gene_name.strip())
          else:
            gene_name="-" #if there is no gene name
          strand=transcript_strand[transcript]
          orf_coordinate=orf_id.split("|")[2].strip()
          orf_coordinates.append(orf_coordinate) #orf genome coordinates
          longest_orf="N"
          if orf_id.split("|")[1] in transcript_longest_orf_map.keys():
            longest_orf="Y"
          
          #uniprot proteins
          if protein_seq in uniprot.keys():
            if uniprot[protein_seq].split("|")[0]=="sp":
              protein_status="reviewed(Swiss-Prot)"
            elif uniprot[protein_seq].split("|")[0]=="tr":
              protein_status="unreviewed(TrEMBL)"
              
            protein_accession=uniprot[protein_seq].split("|")[1]
            longest_orf_anno[protein_accession]=longest_orf
            protein_description=uniprot[protein_seq].split("|")[2]
            
            orfbiotype_transcript_map.setdefault(protein_accession,[]).append(transcript+"|CDS")
            #annotate ORFs
            get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"annotated","CDS","-",longest_orf,protein_status,orf_metadata_annotated_proteins)
          
          #Refseq/Ensembl proteins
          elif protein_seq in refprot.keys():
            protein_status="-"
            protein_accession=refprot[protein_seq]
            longest_orf_anno[protein_accession]=longest_orf
            protein_description="-"
            orfbiotype_transcript_map.setdefault(protein_accession,[]).append(transcript+"|CDS")
            
            #annotate ORFs
            get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"annotated","CDS","-",longest_orf,protein_status,orf_metadata_annotated_proteins)
    #annotated_proteins close
    
    #write output to metadata and fasta file
    canonical_annotated_orf_map={} #k:protein_seq, v:fasta header annotations
    all_annotated_orf_map={} #k:protein_seq, v:fasta header annotations

    for protein_seq,annotations in orf_metadata_annotated_proteins.items():

      for protein_anno in annotations:
        #from metadata

        protein_accession=protein_anno.split("\t")[0] 
        gene_id=protein_anno.split("\t")[1]
        gene_name=protein_anno.split("\t")[2]
        protein_description=protein_anno.split("\t")[3]
        transcript=protein_anno.split("\t")[4]
        orf_coordinate=protein_anno.split("\t")[8]
        longest_orf=protein_anno.split("\t")[14]
       
        if args[6]=="canonical":
         if longest_orf=="Y":
            fw_proteomedb_metadata.write(protein_anno)
            fa_header=protein_accession+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript+"|"+protein_description
            #print(fa_header)
            canonical_annotated_orf_map.setdefault(protein_seq,[]).append(fa_header)
        elif args[6]=="all":
          fw_proteomedb_metadata.write(protein_anno)
          fa_header=protein_accession+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript+"|"+protein_description
          #print(fa_header)
          all_annotated_orf_map.setdefault(protein_seq,[]).append(fa_header)

    if args[6]=="canonical":
      #generate fasta database
      formatFastaheader(canonical_annotated_orf_map,fw_proteomedb)
    elif args[6]=="all":
      #generate fasta database
      formatFastaheader(all_annotated_orf_map,fw_proteomedb)

    #Novel proteins
    AN=Annotations()
    orf_annotation_map={} #stores ORF annotations k:temparory ORF_id, annotations
    counter=1 #counter for temparary ORF ids
    protein_description="-"

    for i in cdhit_clustering_output:
      if '*' in i.strip():
        if ">Bambu" in i.strip() or ">ENS" in i.strip():
          pattern = r'>ENS.*P.*'
          if re.search(pattern, i.strip()): #remove representative sequences of known proteins
            continue
          else:
            longest_seq_orf_id=i.strip().split(">")[1].split("...")[0] #orf_id of Longest representative ORF in cluster
            protein_status="-"
            protein_accession="ORF_"+str(counter)
            protein_seq = unannotated_protein_coordinates[longest_seq_orf_id] #Protein sequence
            #check if protein annotated in openprot
            if protein_seq in openprot.keys():
              openprot_id=openprot[protein_seq]
            else:
              openprot_id="-"
            
            for orf_id in unannotated_proteins[protein_seq]: #accessing all ORFs coordinates for a given protein sequence
              transcript=orf_id.split("|")[0].split("_")[0]
              if transcript in transcript_biotypes.keys():
                transcript_biotype=transcript_biotypes[transcript]
                transcript_coordinates=transcript_genome_coordinates[transcript]
                gene_id=transcript_gene_id_map[transcript]

                if re.match(r'.*\S.*',transcript_gene_name_map[transcript].strip()):
                  gene_name=transcript_gene_name_map[transcript].strip()
                else:
                  gene_name="-"
                  
                longest_orf="N"
                if orf_id.split("|")[0] in transcript_longest_orf_map.keys():
                  longest_orf="Y"
    
                strand=transcript_strand[transcript]
                orf_coordinate=orf_id.split("|")[1].strip()
                
                if transcript_biotype =="protein_coding" and transcript in utr_coordinates.keys():
                  
                  utr_orf=AN.UTRAnnotations(transcript,utr_coordinates[transcript],orf_coordinate,strand,cds_coordinates[transcript])
                  if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf: #ORF overlap with UTR region
                    if len(protein_seq) < int(args[7]): #for short uORFs, remove CDS annotation
                      utr_orf=utr_orf.split(":")[0]
                      #annotate ORF
                      get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated",utr_orf,openprot_id,longest_orf,protein_status,orf_annotation_map)
                    else:
                      #annotate ORF
                      get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated",utr_orf,openprot_id,longest_orf,protein_status,orf_annotation_map)
                   
                      
                elif transcript_biotype !="protein_coding":
                  
                  if transcript in utr_coordinates.keys(): #exceptions: some transcripts have UTR but they are not protein coding
                    
                    utr_orf=AN.UTRAnnotations(transcript,utr_coordinates[transcript],orf_coordinate,strand,cds_coordinates[transcript])
                    if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf: #ORF overlap with UTR region
                      if len(protein_seq) < int(args[7]):
                        utr_orf=utr_orf.split(":")[0]
                        #annotate ORF
                        get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated",utr_orf,openprot_id,longest_orf,protein_status,orf_annotation_map)
                      else:
                        #annotate ORF
                        get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated",utr_orf,openprot_id,longest_orf,protein_status,orf_annotation_map)
                        
                    else:#if ORF is not in UTR regions
                      gene_overlap=AN.isIntergenic(orf_coordinate,protein_coding_gene_coordinates)

                      if gene_overlap=="gene_overlap":
                        #annotate ORF
                        get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated","gene_overlap",openprot_id,longest_orf,protein_status,orf_annotation_map)
                      else:
                        #annotate ORF
                        get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated","intergenic",openprot_id,longest_orf,protein_status,orf_annotation_map)
                  
                  #noncoding RNA transcripts. transcript doesn't have UTRs
                  else:
                    gene_overlap=AN.isIntergenic(orf_coordinate,protein_coding_gene_coordinates)
                    if gene_overlap=="gene_overlap":
                      #annotate ORF
                      get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated","gene_overlap",openprot_id,longest_orf,protein_status,orf_annotation_map)
                    else:
                      #annotate ORF
                      get_protein_annotation(transcript,gene_id,gene_name,protein_description,var_transcript_ORF_map,protein_seq,protein_accession,strand,transcript_biotype,transcript_coordinates,orf_coordinate,"unannotated","intergenic",openprot_id,longest_orf,protein_status,orf_annotation_map)

            #for loop closed
            counter=counter+1

    cdhit_clustering_output.close()
    
    #write output of novel ORFs
    canonical_novel_orf_ids=[]
    all_novel_orf_ids=[]
    
    #store temparary ids of canonical or all novel ORFs according to options provided
    for protein_seq,annotations in orf_annotation_map.items():
      for protein_annotation in annotations:
        longest_orf=protein_annotation.split("\t")[12]
        accession=protein_annotation.split("\t")[0].replace("_var","")
        if args[6]=="canonical":
          if longest_orf=="Y":
            if accession not in canonical_novel_orf_ids:
              canonical_novel_orf_ids.append(accession)
        if args[6]=="all":
          if accession not in all_novel_orf_ids:
            all_novel_orf_ids.append(accession)
    #loop closed

    #replace the temporary ORFs ids based on their index
    canonical_novel_orf_map={} #k:protein_seq, v:fasta header annotations
    all_novel_orf_map={} #k:protein_seq, v:fasta header annotations
    
    for protein_seq,annotations in orf_annotation_map.items():
      for protein_annotation in annotations:
        accession=protein_annotation.split("\t")[0].replace("_var","")
        gene_id=protein_annotation.split("\t")[1]
        gene_name=protein_annotation.split("\t")[2]
        protein_description=protein_annotation.split("\t")[3]
        transcript=protein_annotation.split("\t")[4]
        orf_coordinate=protein_annotation.split("\t")[8]
        longest_orf=protein_annotation.split("\t")[14]
        
        if args[6]=="canonical":
          if longest_orf=="Y":
            new_accesion="ORF_"+str(canonical_novel_orf_ids.index(accession)+1) # +1 as list index starts with zero
            fa_header=new_accesion+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript
            
            canonical_novel_orf_map.setdefault(protein_seq,[]).append(fa_header)
            revised_annotations=protein_annotation.replace(accession,new_accesion)
            fw_proteomedb_metadata.write(revised_annotations)
            
        if args[6]=="all":
          new_accesion="ORF_"+str(all_novel_orf_ids.index(accession)+1) # +1 as list index starts with zero
          fa_header=new_accesion+"|"+orf_coordinate+"|"+gene_id+"|"+gene_name+"|"+transcript
          all_novel_orf_map.setdefault(protein_seq,[]).append(fa_header)
          revised_annotations=protein_annotation.replace(accession,new_accesion)
          fw_proteomedb_metadata.write(revised_annotations)
    #loop closed
    
    #write fasta output
    if args[6]=="canonical":
      #generate fasta database
      formatFastaheader(canonical_novel_orf_map,fw_proteomedb)
    elif args[6]=="all":
      #generate fasta database
      formatFastaheader(all_novel_orf_map,fw_proteomedb)
      
    #remove intermidiate files
    if os.path.isfile(outdir+"orf_temp.txt"):
      subprocess.run(['rm', '-rf', outdir+"orf_temp.txt"], check=True)
    
    if os.path.isfile(outdir+"cdhit_out"):
      subprocess.run(['rm', '-rf', outdir+"cdhit_out"], check=True)
    subprocess.run(['rm', '-rf', outdir+"cdhit_out.clstr"], check=True)

if __name__ == "__main__":
    main()
