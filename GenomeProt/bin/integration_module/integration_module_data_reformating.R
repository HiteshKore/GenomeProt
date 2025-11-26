library(dplyr)
library(stringr)
library(tidyr)
library(optparse)
library(readr)

#usage: Rscript integration_module_data_reformating.R -d <peptide files directory> -s <Proteomics search tool:Spectronaut,FragPipe,FragPipe_quant

#peptide redundency removal
reformat_spectronaut_data<-function(peptide_file,Dataset_id){
  #peptide_file="Combined_peptide_data.tsv"
  peptide_data=read_tsv(peptide_file)
  #consider columns with MS2 intensities and Qvalue
  peptide_data_df<-peptide_data%>%dplyr::select(contains("PG."),contains("PEP.AllOccurringProteinAccessions"),contains("EG.PrecursorId")) %>%
  #replace NA with dummy number 999 in columns with Qvalue
   mutate(Peptide = sapply(str_split(str_replace_all(EG.PrecursorId,"\\[.*?\\]", ""), "_"), `[`, 2))
  
  metadata<-peptide_data_df%>%dplyr::select(Peptide,PEP.AllOccurringProteinAccessions)%>%unique()%>%mutate(Dataset_id=Dataset_id)
  

  return(metadata)
  
}

reformat_fragpipe_data<-function(peptide_file,Dataset_id){
  
  peptide_data=read_tsv(peptide_file)
  peptide_data_flt<-peptide_data%>%dplyr::select(Peptide,Protein,`Mapped Proteins`)
  
  peptide_data_flt <- peptide_data_flt %>%
    mutate(
      # Extract accession from Protein (middle piece between | |)
      Protein = str_extract(Protein, "(?<=\\|)[^|]+(?=\\|)"),
      
      # For Mapped Proteins:
      `Mapped Proteins` = str_split(`Mapped Proteins`, ",\\s*"),
      `Mapped Proteins` = lapply(`Mapped Proteins`, function(x){
        if(all(is.na(x))) return("")
        str_extract(x, "(?<=\\|)[^|]+(?=\\|)")
      }),
      
      # Collapse vector back into comma-separated string
      `Mapped Proteins` = sapply(`Mapped Proteins`, function(x){
        if(all(is.na(x))) return("")
        paste(x, collapse = ",")
      })
    )%>%rename(Mapped_proteins=`Mapped Proteins`)%>%mutate(Dataset_id=Dataset_id)
return(peptide_data_flt)
}



#reformat fragpipe quant data

reformat_fragpipe_quant_data<-function(peptide_file,Dataset_id){
  
  peptide_data=read_tsv(peptide_file)
  
  peptide_data_flt<-peptide_data%>%dplyr::select(`Stripped.Sequence`,Protein.Ids,`All Mapped Proteins`)
  colnames(peptide_data_flt)<-c("Peptide","Protein","Mapped_proteins")
  peptide_data_flt<-peptide_data_flt%>%mutate(
    # Extract accession from Protein (middle piece between | |)
    Protein = str_extract(Protein, "(?<=\\|)[^|]+(?=\\|)"),
    
    # For Mapped Proteins:
    Mapped_proteins= str_split(Mapped_proteins, ","),
    Mapped_proteins = lapply(Mapped_proteins, function(x){
      if(all(is.na(x))) return("")
      str_extract(x, "(?<=\\|)[^|]+(?=\\|)")
    }),
    
    Mapped_proteins = mapply(function(mapped, prot){
      if(is.null(mapped) || length(mapped) == 0) return("")
      
      mapped <- mapped[!is.na(mapped)]     # remove NAs
      mapped <- trimws(mapped)             # remove spaces
      mapped <- mapped[mapped != prot]     # remove the protein accession
      
      if(length(mapped) == 0) "" else paste(mapped, collapse = ",")
    }, Mapped_proteins, Protein),
    
    # Ensure character column
    Mapped_proteins = as.character(Mapped_proteins),
    Dataset_id=Dataset_id
   )
  
  return(peptide_data_flt)
  
}


# define options
option_list = list(
  make_option(c("-d", "--data_directory"), type="character", default=NULL,
              help="peptide data directory", metavar="character"),
  make_option(c("-s", "--search_tool"), type="character", default=NULL,
              help="search tool used to generate proteomic data: Spectronaut,FragPipe,FragPipe_quant", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

peptide_data_dir=opt$data_directory
data_input=opt$search_tool

peptide_files=list.files(peptide_data_dir)


peptide_results<-list()

for (fn in peptide_files) {
  
    file=paste0(peptide_data_dir,fn)
    if(file.exists(file)){
      
      data_type=str_replace(fn,".tsv","")
      if(data_input =="Spectronaut"){
        results<-reformat_spectronaut_data(peptide_file = file,Dataset_id=data_type)
      } else if(data_input =="FragPipe"){
        
        results<-reformat_fragpipe_data(file,data_type)
      }else if(data_input =="FragPipe_quant"){
      
        results<-reformat_fragpipe_quant_data(file,Dataset_id=data_type)
    }
    peptide_results[[data_type]] <- results
    #peptide_results_quant_flt[[data_type]] <- results$peptide_quant_flt
}

}



combined_peptide_results <- do.call(rbind, peptide_results)

if(data_input =="Spectronaut"){
  
  combined_peptide_results_uniq<-combined_peptide_results %>%
    separate_rows(PEP.AllOccurringProteinAccessions, sep = ";") %>%
    group_by(Peptide) %>%
    summarise(
      PEP.AllOccurringProteinAccessions = paste(unique(trimws(PEP.AllOccurringProteinAccessions)), collapse = ","),
      Dataset_id = paste(unique(trimws(Dataset_id)), collapse = ","),
      .groups = 'drop'
    )
  
  uniq_peptides<- combined_peptide_results_uniq%>%mutate(Protein = sapply(str_split(PEP.AllOccurringProteinAccessions, ","), `[`, 1)) %>%
    dplyr::rename(Mapped_proteins=PEP.AllOccurringProteinAccessions)%>%
    rowwise() %>%
    mutate(Mapped_proteins = paste(str_split(Mapped_proteins,",")[[1]][-1], collapse = ","))%>%ungroup%>%
    dplyr::select("Peptide","Protein","Mapped_proteins","Dataset_id")%>%na.omit()
  write_tsv(uniq_peptides,paste0(peptide_data_dir,"peptide_data.tsv"))
  
}else if(data_input =="FragPipe" ||  data_input =="FragPipe_quant" ){
  
  combined_peptide_results_uniq<-combined_peptide_results %>%
    separate_rows(Mapped_proteins, sep = ",") %>%
    group_by(Peptide) %>%
    summarise(
      Protein =paste(unique(trimws(Protein)), collapse = ","),
      Mapped_proteins = paste(unique(trimws(Mapped_proteins)), collapse = ","),
      Dataset_id = paste(unique(trimws(Dataset_id)), collapse = ","),
      .groups = 'drop'
    )
  write_tsv(combined_peptide_results_uniq,paste0(peptide_data_dir,"peptide_data.tsv"))
}







