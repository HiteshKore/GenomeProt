import pysam
import sys
import os
import re
import vcfpy
from vcfpy import Substitution
from collections import Counter
#This script parses the VCF file (v4.2) and segregates homozygous and heterozygous variants into two files.
#It can also handle a merged VCF file generated using bcftools merge or gatk joint variant calling. It simplifies the VCF by
#taking allele with maximum depth in case of mutiple variant alleles. It modifies the allele attributes accordingly  
#Usage: python3 vcfparser.py <vcf file> <outdir>
args=sys.argv


# Convert the record to a plain VCF line
def record_to_vcf_line(record): #returns atributes for first sample in case of multiple samples
    # Extract basic fields
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF

    alts = ",".join([alt.value for alt in record.ALT])  # Extract alternative alleles
    
    qual = record.QUAL if record.QUAL is not None else '.'  # Handle missing quality
    
    filters = ",".join(record.FILTER) if record.FILTER else '.'
    
    # Format INFO field (empty in this case)
    info = '.'
    
    # Format FORMAT field
    format_field = ":".join(record.FORMAT)
    
    # Process sample calls
    call_fields = []
    for call in record.calls:
      call_data = [call.data.get(field, '.') for field in record.FORMAT] #fetch the values for GT:AD:DP:GQ:PL and store in list. If value does not exist it will return '.'
      call_fields.append(":".join(map(str, call_data)))
    
    # Join all parts to form the VCF line
    vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alts}\t{qual}\t{filters}\t{info}\t{format_field}\t{call_fields[0]}" #Just to simplify vcf, it will only consider the attributes from first sample.
    #It doesn't affect the downstream calculations.
    
    vcf_line=vcf_line.replace("[","").replace("]", "").replace(" ", "")
   
    return vcf_line
    
#get variant with maximum depth in case of multiple variants at genomic locus
def getMostFrequentVariant(record,likely_genotype,total_dp,total_ref_ale_dp):
  alt_ale_cnt=len(record.ALT)
  vcf_record=""
  alt_alleles = [alt.value for alt in record.ALT]
  
  varmap={} #k:variant, v:total depth 
  snv_vars={} #stores SNVs. k: variant, total depth
  
  for var in record.ALT:
    #if var.value !="*": 
    varmap[var.value]=0 #initiate value for each alt allele
      
  for call in record.calls:
    # Extract AD value
    ad = call.data.get("AD",0)
    
    if not isinstance(ad, int) and len(ad)>0 and ad is not None: #handles format change in joint vcf produced by GATk and bcftool merge
      ad.pop(0)
      #loop through alternate allele depths
      for ind in range(0,len(ad)):
        alt_ale=list(varmap.keys())[ind]
        
        if ad[ind] is not None:
          varmap[alt_ale]=varmap[alt_ale]+int(ad[ind]) #calculate depth for each allele
          
  
  #for loop close
  for k in varmap.keys():
    #only store SNVs
    if len(k)==1 and k!="*":
      snv_vars[k]=varmap[k]
  
  if len(snv_vars.keys())>0: #if there is at least one SNV
    alt_alele = max(snv_vars, key=snv_vars.get) #variant with maximum depth
    ale_alele_dp=max(list(snv_vars.values())) #alt allele depth
    
    alt_ale_dp_nw=[total_ref_ale_dp,ale_alele_dp] #combine reference allele depth along with alt allele depth
    # Replacing the ALT list
    record.ALT = [Substitution(type_='SNV', value=alt_alele)]
    
    for call in record.calls:
      # Replace AD with the new value
      call.data['AD'] = alt_ale_dp_nw #update alternate depth
      call.data['DP'] = total_dp #update total depth
      call.data['GT'] = likely_genotype #update overall genotype
      
    vcf_record=record_to_vcf_line(record)
    

  return vcf_record
  

#####################################################

vcf_reader = vcfpy.Reader.from_path(args[1]) #vcf file input to vcfpy

file_name=os.path.basename(args[1]) #input file name

#output files
fw=open( args[2]+file_name.replace(".vcf","_homozygous.vcf"),"w")
fw1=open(args[2]+file_name.replace(".vcf","_heterozygous.vcf"),"w")

vcf_fh=open(args[1],"r")
for rec in vcf_fh:
  if rec.strip().startswith('#'):
    if rec.strip().startswith('#CHROM'):
      columns=rec.strip().split("\t")
      format_index=columns.index("FORMAT")
      subset_list = columns[0:format_index +2]
      modified_header=("\t".join(subset_list))
      fw.write(modified_header+"\n")
      fw1.write(modified_header+"\n")
    else:
      fw.write(rec.strip()+"\n")
      fw1.write(rec.strip()+"\n")
  else:
    break
  
vcf_fh.close()

#read vcf_reader
for record in vcf_reader:
  if len(record.REF)==1: #check reference allele is SNV
    depths=[]
    gts=[]
    samples=[] #samples in the analysis
    reference_allele_dp=[]
   
    
    for call in record.calls:
      dp=call.data.get('DP', 0)
      samples.append(call.sample)
      #depth
      if dp is not None:#consider depths with ale allele predicted
        depths.append(dp)
      #genotype
      if '.' not in call.data.get('GT', 0):
        gts.append(call.data.get('GT', 0)) #append genotype to gts 
        
    
      if not isinstance(call.data.get('AD', 0), int) and len(call.data.get('AD',0))>0 and call.data.get('AD',0)[0] is not None: #handles format change in joint vcf produced by GATk and bcftool merge
        #reference allele depth
        reference_allele_dp.append(int(call.data.get('AD',0)[0]))
    
    #for loop closed
    
    #calculate fraction of samples genotype detected in 
    per_samples_with_genotype=len(gts)/len(samples)*100
    
    
    if per_samples_with_genotype>=60: #if genotype is detected in at least 60% samples
      DP=sum(depths)
      unique_gts_map = dict(Counter(gts))
      
      total_ref_allele_dp=sum(reference_allele_dp) #overall depth of reference allele
      
      most_likely_genotype= max(unique_gts_map, key=unique_gts_map.get) #most likely genotype based on majority rule
      
      if DP >=8 and int(record.QUAL)>10:
        if int(re.split(r"[\/|]", most_likely_genotype)[0])>0: # homozygous records '|' phased genotype: when maternal and paternal allele known
          if total_ref_allele_dp <10: #there might be some misalignments if reference alele is less than 10, mutation still will be considered homozygous
            
            vcf_rec=getMostFrequentVariant(record,most_likely_genotype,DP,total_ref_allele_dp)
            if vcf_rec.strip(): #remove empty lines
              fw.write(vcf_rec+"\n")
        else:
          if total_ref_allele_dp > 0: #check if reference allele depth is not zero
            
            vcf_rec=getMostFrequentVariant(record,most_likely_genotype,DP,total_ref_allele_dp)
            if vcf_rec.strip(): #remove empty lines
              fw1.write(vcf_rec+"\n")
