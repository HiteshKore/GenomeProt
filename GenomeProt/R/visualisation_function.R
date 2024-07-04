
# define visualisation function
plot_gene <- function(gene_name, tx_res, pep_res, orf_res, txcounts=NA, pepcounts=NA, min_intron_len=500) {
  
  # test vars
  # gene_name <- "ENSG00000169220.18"
  # tx_res <- res_tx_import
  # pep_res <- res_pep_import
  # orf_res <- res_ORF_import
  # min_intron_len=500
  # txcounts <- countstm
  # pepcounts <- countspm
  
  if (!(gene_name %in% tx_res$gene_id)) {
    print("Provided gene was not found in data")
    stop()
  }
  
    # filter for selected gene
    tx_res <- tx_res %>% 
      dplyr::filter(!is.na(gene_id), gene_id == gene_name) %>% 
      mutate(feature_type = "Transcripts",
             peptide_type = "Transcripts",
             ORF_id = NA) %>% 
      dplyr::select(seqnames,start,end,strand,type,gene_id,transcript_id,feature_type,peptide_type,exon_number,ORF_id)
    
    # filter for selected gene
    pep_res$transcript_id <- pep_res$peptide
    pep_res <- pep_res %>% 
      dplyr::filter(!is.na(gene_id), gene_id == gene_name) %>% 
      separate(naming, into="tx_id", sep="_")
    
    pep_res <- pep_res %>% 
      mutate(feature_type = "Peptides",
             peptide_type = pep_map_status) %>% 
      group_by(transcript_id) %>% 
      mutate(ORF_id = case_when(
        exon_number == 1 ~ PID,
        TRUE ~ NA)) %>% 
      ungroup() %>% 
      separate(ORF_id, into="ORF_id", sep="\\|") %>% 
      dplyr::select(seqnames,start,end,strand,type,gene_id,transcript_id,tx_id,feature_type,peptide_type,exon_number,ORF_id)
    
    # filter for selected gene
    orf_res$transcript_id <- orf_res$tx_id
    orf_res$ORF_id <- NULL
    
    orf_res <- orf_res %>% 
      dplyr::filter(tx_id %in% pep_res$tx_id) %>% 
      mutate(type = "CDS",
             feature_type = "Transcripts",
             gene_id = gene_name,
             peptide_type = "Transcripts") %>% 
      group_by(transcript_id) %>% 
      mutate(ORF_id = case_when(
        exon_number == 1 ~ PID,
        TRUE ~ NA)) %>% 
      ungroup() %>% 
      separate(ORF_id, into="ORF_id", sep="\\|") %>% 
      dplyr::select(seqnames,start,end,strand,type,gene_id,transcript_id,feature_type,peptide_type,exon_number,ORF_id)
    
    # filter for exons only in transcripts gtf
    gtf_exons <- tx_res %>% 
      dplyr::filter(type == "exon" & transcript_id %in% pep_res$tx_id)
    
    pep_res$tx_id <- NULL
    
    # combine peptides, transcripts and ORFs
    gtf_to_plot <- rbind(pep_res, orf_res, gtf_exons)
    
    # factor and set levels
    gtf_to_plot$feature_type <- factor(gtf_to_plot$feature_type, levels=c('Peptides', 'Transcripts'))
    gtf_to_plot$peptide_type <- factor(gtf_to_plot$peptide_type, levels=c('low', 'medium', 'high', 'Transcripts'))
    
    # add a '*' label when a peptide is uniquely mapped
    gtf_to_plot <- gtf_to_plot %>% 
      arrange(peptide_type) %>% 
      mutate(ORF_id = case_when(
        feature_type == "Transcripts" ~ ORF_id,
        feature_type == "Peptides" & peptide_type == "high" & exon_number == 1 ~ ORF_id,
        feature_type == "Peptides" & peptide_type != "high" & exon_number != 1 ~ NA
      ))
    
    pep_gtf_to_plot <- gtf_to_plot %>% 
      dplyr::filter(feature_type == "Peptides")
    
    tx_gtf_to_plot <- gtf_to_plot %>% 
      dplyr::filter(feature_type == "Transcripts")
    
    # define CDS as the ORFs
    gtf_cds <- tx_gtf_to_plot %>% 
      dplyr::filter(type == "CDS") %>% 
      mutate(ORF_id = NA)
    
    # number of unique transcripts and peptides for panel heights
    n_tx <- length(unique(tx_res$transcript_id))
    n_pep <- length(unique(pep_res$transcript_id))
    
    # before plotting, check for consistent transcripts and peptides
    if (missing(txcounts) & missing(pepcounts)) {
      
      #skip
      
    } else {
      
      # add quantitative heatmap info
      # filter transcript counts for gene of interest
      txcounts <- txcounts %>%
        dplyr::filter(transcript_id %in% tx_gtf_to_plot$transcript_id)
      # filter peptide counts for gene of interest
      pepcounts <- pepcounts %>%
        dplyr::filter(peptide %in% pep_gtf_to_plot$transcript_id)
      
      txcounts$sample_id <- as.factor(as.character(txcounts$sample_id))
      pepcounts$sample_id <- as.factor(as.character(pepcounts$sample_id))
      
      tx_gtf_to_plot <- tx_gtf_to_plot %>% 
        dplyr::filter(transcript_id %in% txcounts$transcript_id)
      pep_gtf_to_plot <- pep_gtf_to_plot %>% 
        dplyr::filter(transcript_id %in% pepcounts$peptide)
      
    }
    
    # plot peptide and transcript tracks
    gtf_tx_output <- tx_gtf_to_plot %>% 
      ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
      geom_intron(data = to_intron(tx_gtf_to_plot, "transcript_id"),
                  aes(strand = strand), arrow.min.intron.length = min_intron_len) +
      geom_range(aes(fill = peptide_type), show.legend = FALSE, height = 0.5) +
      geom_range(data = gtf_cds, height = 0.75, fill = "#295D9B") +
      ylab("") +
      xlab(paste0(unique(tx_gtf_to_plot$seqnames), " ", unique(tx_gtf_to_plot$gene_id))) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text.y = element_blank()) +
      scale_fill_manual(values = c("low" = "#D3D3D3", "medium" = "#818589", "high" = "orangered2", "Transcripts" = "#9FC9FB")) +
      geom_text_repel(aes(x = start, label = ORF_id), size = 3, nudge_y = 0.5, min.segment.length = Inf)
    
    xlimits <- c(layer_scales(gtf_tx_output)$x$range$range)
    
    gtf_pep_output <- pep_gtf_to_plot %>% 
      ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
      geom_intron(data = to_intron(pep_gtf_to_plot, "transcript_id"),
                  aes(strand = strand), arrow.min.intron.length = min_intron_len) +
      geom_range(aes(fill = peptide_type), show.legend = FALSE, height = 0.5) +
      ylab("") +
      xlab("") +
      xlim(xlimits) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_fill_manual(values = c("low" = "#D3D3D3", "medium" = "grey16", "high" = "orangered2", "Transcripts" = "#9FC9FB")) +
      geom_text_repel(aes(x = start, label = ORF_id), size = 3, nudge_y = 0.5, min.segment.length = Inf)
    
    # return peptide and transcript tracks if no quant data is provided
    if (missing(txcounts) & missing(pepcounts)) {
      
      pep_vis_plot <- gtf_pep_output + gtf_tx_output + 
        plot_layout(nrow = 2, ncol = 1, heights = c(n_pep, n_tx))
      
    } else {
      
      
      # plot peptide heatmap
      pepcounts_output <- pepcounts %>%
        ggplot(aes(x = sample_id, y = peptide, fill = count)) +
        geom_tile() +
        scale_fill_viridis_b(n.breaks = 8, option = "D") +
        labs(x = "", y = "", fill = "VSN Intensity") +
        theme_bw() +
        theme(strip.background = element_rect(color = "black", fill = "white", size = 0, linetype = "solid"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              #axis.text.x = element_blank(),
              #axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.text = element_text(size=6))
      
      # plot transcript heatmap
      txcounts_output <- txcounts %>%
        ggplot(aes(x = sample_id, y = transcript_id, fill = log2(count+1))) +
        geom_tile() +
        scale_fill_viridis_b(n.breaks = 8, option = "D") +
        labs(x = "Sample", y = "", fill = "log2(TPM+1)") +
        theme_bw() +
        theme(strip.background = element_rect(color = "black", fill = "white", size = 0, linetype = "solid"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.text = element_text(size=6))
      
      pep_vis_plot <- gtf_pep_output + pepcounts_output + gtf_tx_output + txcounts_output + 
        plot_layout(nrow = 2, ncol = 2, widths = c(1.5, 1), heights = c(n_pep, n_tx))
      
    }
  return(pep_vis_plot)
}

