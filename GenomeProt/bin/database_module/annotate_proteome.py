#This script removes the redundant ORFs and considers longest ORF if ORF is part of longer ORF.
#It annotates them based on their location on the genome/transcript
#cd-hit (First install cd-hit commond line tool in linux (command to install:'conda install bioconda/label/cf201901::cd-hit'). Next, install py-cdhit library using  'pip install py-cdhit' command)
#Usage:python3 annotate_proteome.py gencode.vM33.chr_patch_hapl_scaff.annotation_chrX.gtf openprot_uniprotDb_mm.txt  ORFome_aa.txt proteome_database_transcripts.gtf <outdir> <canonical/all> <orf_length> <mutant_protein_db>
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
  if len(args) != 9:
        print("Usage: python3 annotate_proteome.py <reference_gtf> <custom_openprot+uniprot_db> <ORFome_aa.txt> <ORFome_transcripts.gtf> <outdir> <canonical/all> <orf_length> <mutant_protein_db>")
        sys.exit(1)
  else:
    
    #custom openprot+ uniprot annotation database
    RF=refDb()
    uniprot={} #comprise UniProt annotations: k:seq v: trEMBL/reviewed|protein_id| gene description
    openprot={} #comprise OpenProt annotations: k:seq v:protein_id
    refprot={} #comprise reference proteins annotated in uniprot/refseq/ensembl. k:seq v:protein id
    annotated_protein_ids={}
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
    
    #ORFomedb
    unannotated_proteins={} #stores unannotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    unannotated_protein_coordinates={} #stores unannotated protein coordinates as key: k: transcript|orf_number|genomic_coordinates,v:protein sequence
    orfome=open(args[3]) 
    transcript_orf_map={} #stores all orfs for each transcript: k;transcript_id,v:protein_seq|orf_id
    
    annotated_proteins={} #store annotated proteins: k:protein sequence, v:transcript|orf_number|genomic_coordinates
    protein_seq_properties={}#k:protein_sequence,v:physicochemical_properties 
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
            
            
            
          
        #calculate physico-chemical propertise
        SP=SequenceProperties()
        Mol_wt=SP.calculateMolWt(protein_seq)
        IsoPt=SP.calculateIsoElectricPoint(protein_seq)
        HydroP_ind=SP.calculateHydrophobicity(protein_seq)
        Aliphatic_ind=SP.calculateAliphatic_index(protein_seq)
        seq_properties=str(Mol_wt)+"\t"+str(IsoPt)+"\t"+str(HydroP_ind)+"\t"+str(Aliphatic_ind)
        protein_seq_properties[protein_seq]=seq_properties
        
        
        
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
    fw_proteomedb_metadata.write("accession\tgene\tgene_symbol\ttranscript\tstrand\ttranscript_biotype\ttranscript_coordinates\torf_genomic_coordinates\torf_type\tlocalisation\topenprot_id\tprotein_sequence\tlongest_orf_in_transcript\tuniprot_status\tamino_acid_change\tmolecular_weight(kDA)\tisoelectric_point\thydrophobicity\taliphatic_index\n")
        

    gene_ORF_map={} #k: gene_id v: ORF_ids
    longest_orf_anno={} #k:accession v-:Y, N
    orfome_transcripts=[]
    orfbiotype_transcript_map={} #k:orf, v:transcript|biotype
    #annotations of known proteins
    
    for k in annotated_proteins.keys():
      protein_seq=k
      transcripts=[]
      gene_accession=[]
      gene_symbol=[]
      orf_coordinates=[]
      #longest ORF
      transcripts_l=[]
      gene_accession_l=[]
      gene_symbol_l=[]
      orf_coordinates_l=[]

      longest_orf_seq=""
      longest_orf_accession=""
      for orf_id in  annotated_proteins[k]:
        transcript=orf_id.split("|")[1].split("_")[0]
        orfome_transcripts.append(transcript) #transcripts contributing to ORFome
        
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
            
            protein_description=uniprot[protein_seq].split("|")[2]
            
            longest_orf_anno[protein_accession]=longest_orf
            
            orfbiotype_transcript_map.setdefault(protein_accession,[]).append(transcript+"|CDS")
            

            if args[6]=="canonical":
              if longest_orf=="Y":
                longest_orf_accession=protein_accession
                #print(longest_orf_accession)
                orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tannotated\tCDS\t"+"-\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                longest_orf_seq=protein_seq
                
                #longest ORF
                gene_accession_l.append(gene_id)
                gene_symbol_l.append(gene_name)
                transcripts_l.append(transcript)
                orf_coordinates_l.append(orf_coordinate)
                fw_proteomedb_metadata.write(orf_annotation)
                gene_ORF_map.setdefault(gene_id,[]).append(protein_accession+"|"+protein_seq) #store protein accessions for gene_id
                

                
            elif args[6]=="all":
              orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tannotated\tCDS\t"+"-\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
              gene_ORF_map.setdefault(gene_id,[]).append(protein_accession+"|"+protein_seq) #store protein accessions for gene_id
              fw_proteomedb_metadata.write(orf_annotation)
              
          #Refseq/Ensembl proteins
          elif protein_seq in refprot.keys():
            protein_status="-"
            protein_accession=refprot[protein_seq]
            longest_orf_anno[protein_accession]=longest_orf
            orfbiotype_transcript_map.setdefault(protein_accession,[]).append(transcript+"|CDS")
            
            if args[6]=="canonical":
              
              if longest_orf=="Y":
                longest_orf_accession=protein_accession
                orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tannotated\tCDS\t"+"-\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                longest_orf_seq=protein_seq
                #longest ORF
                gene_accession_l.append(gene_id)
                gene_symbol_l.append(gene_name)
                transcripts_l.append(transcript)
                orf_coordinates_l.append(orf_coordinate) 
                fw_proteomedb_metadata.write(orf_annotation)
                gene_ORF_map.setdefault(gene_id,[]).append(protein_accession+"|"+protein_seq) #store protein accessions for gene_id
                
            elif args[6]=="all":
              orf_annotation=protein_accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tannotated\tCDS\t"+"-\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
              fw_proteomedb_metadata.write(orf_annotation)
              gene_ORF_map.setdefault(gene_id,[]).append(protein_accession+"|"+protein_seq) #store protein accessions for gene_id
              
  
      #fasta output
      orf_coordinates=list(set(orf_coordinates))
      transcripts=list(set(transcripts))
      gene_accession=list(set(gene_accession))
      gene_symbol=list(set(gene_symbol))
    
      orf_transcripts_s=",".join(transcripts)
      orf_cooordinates_s=",".join(orf_coordinates)
      orf_genes_accession_s=",".join(gene_accession)
      
      if len(gene_symbol)==0:
        gene_symbol.append("-")
      orf_genes_symbol_s=",".join(gene_symbol)
      
      #longest ORF

      orf_transcripts_l_s=",".join(transcripts_l)
      orf_cooordinates_l_s=",".join(orf_coordinates_l)
      orf_genes_accession_l_s=",".join(gene_accession_l)
      
      if len(gene_symbol_l)>=2 and "-" in gene_symbol_l:
        gene_symbol_l.remove('-')
      orf_genes_symbol_l_s=",".join(gene_symbol_l)

      if args[6]=="canonical":
        if longest_orf_accession!="":
          fw_proteomedb.write(">"+longest_orf_accession+"|CO="+orf_cooordinates_l_s+" GA="+orf_genes_accession_l_s+" GN="+orf_genes_symbol_l_s+" TA="+orf_transcripts_l_s+"\n"+protein_seq+"\n")
      
      elif args[6]=="all":
        fw_proteomedb.write(">"+protein_accession+"|CO="+orf_cooordinates_s+" GA="+orf_genes_accession_s+" GN="+orf_genes_symbol_s+" TA="+orf_transcripts_s+"\n"+protein_seq+"\n")
      

    #Novel proteins
    AN=Annotations()
    orf_annotation_map={} #stores ORF annotations k:temparory ORF_id, annotations
    
    
    counter=1 #counter for temparary ORF ids
  
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
              orfome_transcripts.append(transcript) #transcripts contributing to ORFome
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
                    if len(protein_seq) < int(args[7]):
                      utr_orf=utr_orf.split(":")[0]
                      orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\t"+utr_orf+"\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                      orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                    else:
                      orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\t"+utr_orf+"\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                      orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
    
                elif transcript_biotype !="protein_coding":
                  if transcript in utr_coordinates.keys():
                    utr_orf=AN.UTRAnnotations(transcript,utr_coordinates[transcript],orf_coordinate,strand,cds_coordinates[transcript])
                    if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf: #ORF overlap with UTR region
                      if len(protein_seq) < int(args[7]):
                        utr_orf=utr_orf.split(":")[0]
                        orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\t"+utr_orf+"\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                        orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                      else:
                        orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\t"+utr_orf+"\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                        orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                    else:
                      gene_overlap=AN.isIntergenic(orf_coordinate,protein_coding_gene_coordinates)
                      if gene_overlap=="gene_overlap":
                        orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\tgene_overlap\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                        orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                      else:
                        orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\tintergenic\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                        orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                  #if transcript doesn't have UTR
                  else:
                    gene_overlap=AN.isIntergenic(orf_coordinate,protein_coding_gene_coordinates)
                    if gene_overlap=="gene_overlap":
                      orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\tgene_overlap\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                      orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
                    else:
                      orf_annotation=gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tunannotated\tintergenic\t"+openprot_id+"\t"+protein_seq+"\t"+longest_orf+"\t"+protein_status+"\t-\t"+protein_seq_properties[protein_seq]+"\n"
                      orf_annotation_map.setdefault(protein_accession,[]).append(orf_annotation)
    
            counter=counter+1 # for loop closed


    cdhit_clustering_output.close()
    
    orfome_transcripts=list(set(orfome_transcripts)) #uniq transcripts

        

  ##############
    counter1=1 #counter for all protein isoforms
    counter2=1 #counter for canonical proteins
    counter3=1

    for orf_id in orf_annotation_map.keys():
      protein_accession="ORF_"+str(counter1)
      protein_accession_logest_orf="ORF_"+str(counter2)
      transcripts=[]
      gene_accession=[]
      gene_symbol=[]
      orf_coordinates=[]
      
      #longest ORF
      transcripts_l=[]
      gene_accession_l=[]
      gene_symbol_l=[]
      orf_coordinates_l=[]
      longest_orf_seq=""
      #writing output in metadata file
      for orf_annotations in orf_annotation_map[orf_id]:
        protein_seq=orf_annotations.split("\t")[10].strip()
        #####
        gene_id=orf_annotations.split("\t")[0]
        gene_accession.append(gene_id)
        gene_name=orf_annotations.split("\t")[1]
        gene_symbol.append(gene_name)
        transcript=orf_annotations.split("\t")[2]
        transcripts.append(transcript)
        orf_coordinate=orf_annotations.split("\t")[6]
        orf_coordinates.append(orf_coordinate) #orf genome coordinates
        longest_orf=orf_annotations.split("\t")[11]
        
        orf_type=orf_annotations.split("\t")[8]
        
        
        
        if args[6]=="all":
          longest_orf_anno[protein_accession]=longest_orf
          gene_ORF_map.setdefault(gene_id,[]).append(protein_accession+"|"+protein_seq) #store protein accessions for gene_id
          fw_proteomedb_metadata.write(protein_accession+"\t"+orf_annotations)
          orfbiotype_transcript_map.setdefault(protein_accession,[]).append(transcript+"|"+orf_type)
        
        if args[6]=="canonical":
          longest_orf_anno[protein_accession_logest_orf]=longest_orf
          if longest_orf=="Y":
            gene_ORF_map.setdefault(gene_id,[]).append(protein_accession_logest_orf+"|"+protein_seq) #store protein accessions for gene_id
            longest_orf_seq=protein_seq
            fw_proteomedb_metadata.write(protein_accession_logest_orf+"\t"+orf_annotations)
            orfbiotype_transcript_map.setdefault(protein_accession_logest_orf,[]).append(transcript+"|"+orf_type)
            
            
            
            #longest ORF
            gene_accession_l.append(gene_id)
            gene_symbol_l.append(gene_name)
            transcripts_l.append(transcript)
            orf_coordinates_l.append(orf_coordinate) 
            
          #include UTR ORF in addition to canonical
          elif 'UTR' in orf_type and 'CDS:3UTR' not in orf_type :
            gene_ORF_map.setdefault(gene_id,[]).append(protein_accession_logest_orf+"|"+protein_seq) #store protein accessions for gene_id
            orf_annotations_list=orf_annotations.split('\t')
            longest_orf_seq=protein_seq
            fw_proteomedb_metadata.write(protein_accession_logest_orf+"\t"+"\t".join(orf_annotations_list))
            orfbiotype_transcript_map.setdefault(protein_accession_logest_orf,[]).append(transcript+"|"+orf_type)
            
            gene_accession_l.append(gene_id)
            gene_symbol_l.append(gene_name)
            transcripts_l.append(transcript)
            orf_coordinates_l.append(orf_coordinate)
            

        
      orf_coordinates=list(set(orf_coordinates))
      transcripts=list(set(transcripts))
      gene_accession=list(set(gene_accession))
      gene_symbol=list(set(gene_symbol))
    
      orf_transcripts_s=",".join(transcripts)
      orf_cooordinates_s=",".join(orf_coordinates)
      orf_genes_accession_s=",".join(gene_accession)
      
      if len(gene_symbol)>=2 and "-" in gene_symbol:
        gene_symbol.remove('-')
      orf_genes_symbol_s=",".join(gene_symbol)
      #longest ORF
      if len(gene_symbol_l)>=2 and "-" in gene_symbol_l:
        gene_symbol_l.remove('-')
        
      orf_genes_symbol_l_s=",".join(gene_symbol_l)
      orf_transcripts_l_s=",".join(transcripts_l)
      orf_cooordinates_l_s=",".join(orf_coordinates_l)
      orf_genes_accession_l_s=",".join(gene_accession_l)
      
      
      if args[6]=="canonical":
        if re.match(r'.*\S.*',longest_orf_seq.strip()):
          fw_proteomedb.write(">"+protein_accession_logest_orf+"|CO="+orf_cooordinates_l_s+" GA="+orf_genes_accession_l_s+" GN="+orf_genes_symbol_l_s+" TA="+orf_transcripts_l_s+"\n"+protein_seq+"\n")
          
          counter2=counter2+1
      elif args[6]=="all":
        fw_proteomedb.write(">"+protein_accession+"|CO="+orf_cooordinates_s+" GA="+orf_genes_accession_s+" GN="+orf_genes_symbol_s+" TA="+orf_transcripts_s+"\n"+protein_seq+"\n")
      counter1=counter1+1
  
    #Mutant protein annnotations
    if args[8] !="None" and os.path.isfile(args[8]):
      mutant_protein_fh=open(args[8])
      mutant_protein_map={}
      mutant_protein_transcript_id={}
     

      for i in mutant_protein_fh:
        if i.strip().startswith("transcript"):
          continue
        else:
          transcript=i.strip().split("\t")[0]
          gene_id=transcript_gene_id_map[transcript]
          
          #if transcript in orfome_transcripts: #if mutated transcript part of transcritps contributing to wild type ORFome
          mut_sq=i.strip().split("\t")[1]
          orf_coordinate=i.strip().split("\t")[2]
          
          #Adjust orf coordinates
          if strand=="-":
            orf_coordinate=orf_coordinate.split(":")[0]+":"+str(int(orf_coordinate.split(":")[1].split("-")[0])+3)+"-"+orf_coordinate.split("-")[1]
          elif strand=="+":
            orf_coordinate=orf_coordinate.split("-")[0]+"-"+str(int(orf_coordinate.split("-")[1])-3)

          if gene_id in gene_ORF_map.keys():
            for wt_id_seq in gene_ORF_map[gene_id]:
              wt_sq=wt_id_seq.split("|")[1]
              prot_id=wt_id_seq.split("|")[0]
              # percentage similarity and changed residues 
              SS=SequenceSimilarity()
              similarity, aa_change = SS.calculate_similarity(mut_sq, wt_sq)
              similarity=round(similarity,2)
              if similarity >96 and similarity <100:
                mutant_protein_transcript_id.setdefault(mut_sq,[]).append(transcript+"|"+orf_coordinate)
                mutant_protein_map[mut_sq]=prot_id+"|"+wt_sq
      
      mutant_protein_fh.close()
      
      mutant_proteins={}
      
      for mut_sq,orf_ids in mutant_protein_transcript_id.items():
        mutant_protein_gene_id=[]
        mutant_protein_transcript=[]
        mutant_protein_orf_coodinates=[]

        mutant_protein_gene_name=[]
        mut_prot_acc=""
        
        #calculate physico-chemical propertise
        SP1=SequenceProperties()
        Mol_wt=SP1.calculateMolWt(mut_sq)
        IsoPt=SP1.calculateIsoElectricPoint(mut_sq)
        HydroP_ind=SP1.calculateHydrophobicity(mut_sq)
        Aliphatic_ind=SP1.calculateAliphatic_index(mut_sq)
        seq_properties=str(Mol_wt)+"\t"+str(IsoPt)+"\t"+str(HydroP_ind)+"\t"+str(Aliphatic_ind)
        
        
        for orf_id in list(set(orf_ids)):
          transcript=orf_id.split("|")[0]
          orf_coordinate=orf_id.split("|")[1]
          transcript_biotype=transcript_biotypes[transcript]
          transcript_coordinates=transcript_genome_coordinates[transcript]
          gene_id=transcript_gene_id_map[transcript]
          strand=transcript_strand[transcript]
          
          
          
          if re.match(r'.*\S.*',transcript_gene_name_map[transcript].strip()):
            gene_name=transcript_gene_name_map[transcript].strip()
          else:
            gene_name="-"
            
          
          
          wt_prot_id=mutant_protein_map[mut_sq].split("|")[0]
          wt_prot_sq=mutant_protein_map[mut_sq].split("|")[1]
          
          for orfbiotype_transcript in orfbiotype_transcript_map[wt_prot_id]:#check mutant transcript in the wild type transcript map
            wt_transcript=orfbiotype_transcript.split("|")[0]
            wt_orf_localisation=orfbiotype_transcript.split("|")[1]
            
            
            #print(wt_transcript,wt_orf_localisation)
            if "ORF_" not in wt_prot_id:
              accession="m_"+wt_prot_id
              mutant_proteins[mut_sq]=accession
              protein_status="-"
              if wt_prot_sq in uniprot.keys():
                if uniprot[wt_prot_sq].split("|")[0]=="sp":
                  protein_status="reviewed(Swiss-Prot)"
                elif uniprot[wt_prot_sq].split("|")[0]=="tr":
                  protein_status="unreviewed(TrEMBL)"
                  
              if transcript==wt_transcript and wt_orf_localisation=="CDS":
                orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tvariant\tCDS\t-\t"+mut_sq+"\t"+longest_orf_anno[wt_prot_id]+"\t"+protein_status+"\t"+wt_prot_id+":"+aa_change+"\t"+seq_properties+"\n"
                #print(orf_annotation)
                fw_proteomedb_metadata.write(orf_annotation)
                mutant_protein_gene_id.append(gene_id)
                mutant_protein_gene_name.append(gene_name)
                mutant_protein_transcript.append(transcript)
                mutant_protein_orf_coodinates.append(orf_coordinate)
            
            else:
              if mut_sq not in mutant_proteins.keys():
                accession="mORF_"+str(counter3)
                #print(accession,mut_sq)
                mutant_proteins[mut_sq]=accession
                counter3 = counter3 + 1
              
              if "mORF_" in accession:
                #print(accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tMutant\tPH\t-\t"+protein_seq+"\tuniprot_status\t"+seq_properties)
                
                if transcript in utr_coordinates.keys():
                  utr_orf=AN.UTRAnnotations(transcript,utr_coordinates[transcript],orf_coordinate,strand,cds_coordinates[transcript])
                  if "UTR" in utr_orf and "CDS:3UTR" not in utr_orf: #ORF overlap with UTR region
                    utr_orf=utr_orf.split(":")[0]
                    
                    if transcript==wt_transcript and wt_orf_localisation==utr_orf:
                      mutant_protein_gene_id.append(gene_id)
                      mutant_protein_gene_name.append(gene_name)
                      mutant_protein_transcript.append(transcript)
                      mutant_protein_orf_coodinates.append(orf_coordinate)
                      orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tvariant\t"+utr_orf+"\t-\t"+mut_sq+"\t"+longest_orf_anno[wt_prot_id]+"\t-\t"+wt_prot_id+":"+aa_change+"\t"+seq_properties+"\n"
                      fw_proteomedb_metadata.write(orf_annotation)
                  
                  
                  elif "UTR" in utr_orf:
                    if transcript==wt_transcript and wt_orf_localisation==utr_orf:
                      mutant_protein_gene_id.append(gene_id)
                      mutant_protein_gene_name.append(gene_name)
                      mutant_protein_transcript.append(transcript)
                      mutant_protein_orf_coodinates.append(orf_coordinate)
                      orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tvariant\t"+utr_orf+"\t-\t"+mut_sq+"\t"+longest_orf_anno[wt_prot_id]+"\t-\t"+wt_prot_id+":"+aa_change+"\t"+seq_properties+"\n"
                      #print(orf_annotation "no res")
                      fw_proteomedb_metadata.write(orf_annotation)
                  #else:
                    #orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tMutant\t"+utr_orf+"\t-\t"+protein_seq+"\tuniprot_status\t"+seq_properties+"\n"
                else:
                  gene_overlap=AN.isIntergenic(orf_coordinate,protein_coding_gene_coordinates)
                  if gene_overlap=="gene_overlap":
                    if transcript==wt_transcript and wt_orf_localisation=="gene_overlap":
                      mutant_protein_gene_id.append(gene_id)
                      mutant_protein_gene_name.append(gene_name)
                      mutant_protein_transcript.append(transcript)
                      mutant_protein_orf_coodinates.append(orf_coordinate)
                      orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tvariant\tgene_overlap\t-\t"+mut_sq+"\t"+longest_orf_anno[wt_prot_id]+"\t-\t"+wt_prot_id+":"+aa_change+"\t"+seq_properties+"\n"
                      #print(orf_annotation)
                      fw_proteomedb_metadata.write(orf_annotation)
                  else:
                    if transcript==wt_transcript and wt_orf_localisation=="intergenic":
                      mutant_protein_gene_id.append(gene_id)
                      mutant_protein_gene_name.append(gene_name)
                      mutant_protein_transcript.append(transcript)
                      mutant_protein_orf_coodinates.append(orf_coordinate)
                      orf_annotation=accession+"\t"+gene_id+"\t"+gene_name+"\t"+transcript+"\t"+strand+"\t"+transcript_biotype+"\t"+transcript_coordinates+"\t"+orf_coordinate+"\tvariant\tintergenic\t-\t"+mut_sq+"\t"+longest_orf_anno[wt_prot_id] +"\t-\t"+wt_prot_id+":"+aa_change+"\t"+seq_properties+"\n"
                      fw_proteomedb_metadata.write(orf_annotation)
          
          
        #loop close
        if len(list(set(mutant_protein_transcript))) >= 1:
          orf_cooordinates_s=",".join(list(set(mutant_protein_orf_coodinates)))
          orf_genes_accession_s=",".join(list(set(mutant_protein_gene_id)))
          orf_genes_symbol_s=",".join(list(set(mutant_protein_gene_name)))
          orf_transcripts_s=",".join(list(set(mutant_protein_transcript)))
          fw_proteomedb.write(">"+accession+"|CO="+orf_cooordinates_s+" GA="+orf_genes_accession_s+" GN="+orf_genes_symbol_s+" TA="+orf_transcripts_s+"\n"+mut_sq+"\n")


      mutant_protein_fh.close()
     

    #remove intermidiate files
    if os.path.isfile(outdir+"orf_temp.txt"):
      subprocess.run(['rm', '-rf', outdir+"orf_temp.txt"], check=True)
    
    if os.path.isfile(outdir+"cdhit_out"):
      subprocess.run(['rm', '-rf', outdir+"cdhit_out"], check=True)
    subprocess.run(['rm', '-rf', outdir+"cdhit_out.clstr"], check=True)



if __name__ == "__main__":
    main()
