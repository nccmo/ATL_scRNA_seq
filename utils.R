
# RP[S/L], TR[A/B/D/G]V, IG[K/L/H]V exclusion
rps_trv_ig_exclusion <- function(input) {
  tmp <- input[!str_detect(rownames(input), "RPS") & 
                 !str_detect(rownames(input), "RPL") &
                 !str_detect(rownames(input), "TRAV") &
                 !str_detect(rownames(input), "TRBV") &
                 !str_detect(rownames(input), "TRDV") &
                 !str_detect(rownames(input), "TRGV") &
                 !str_detect(rownames(input), "IGKV") &
                 !str_detect(rownames(input), "IGLV") &
                 !str_detect(rownames(input), "IGHV"),]
  return(tmp)
}

# cutoff_filtering function
cutoff_filtering <- function(input, min_pct_cutoff = 0.1, min_logfc = 0.25, p_val = 0.01) {
  tmp <- input[input$pct.1 > min_pct_cutoff &
                 input$avg_log2FC > min_logfc &
                 input$p_val_adj < p_val,]
  return(tmp)
}

rna_deg_add_info <- function(
  deg_list,
  object
) {
  # average expression
  tmp <- AverageExpression(object, assays = "RNA", features = rownames(deg_list))$RNA
  colnames(tmp) <- paste(str_replace_all(colnames(tmp), " ", "_"), "_avg_exp", sep = "")
  
  # pct
  for (i in levels(object)){
    tmp_cells <- Cells(subset(object, ident = i))
    
    tmp[,paste(str_replace_all(i, " ", "_"), "_pct", sep = "")] <- round(
      rowSums(x = object$RNA@data[rownames(deg_list), tmp_cells] > 0, na.rm = TRUE) /
        rowSums(!is.na(object$RNA@data[rownames(deg_list), tmp_cells])),
      digits = 3
    )
  }
  
  deg_list <- deg_list %>% 
    as_tibble(rownames = "Gene")
  tmp <- tmp %>% 
    as_tibble(rownames = "Gene")
  
  output_df <- deg_list %>% 
    left_join2(tmp)
  return(output_df)
}

adt_deg_add_info <- function(
  deg_list,
  object
) {
  # average expression
  tmp <- AverageExpression_na_removal(object, assays = "ADT", features = rownames(deg_list))$ADT
  colnames(tmp) <- paste(str_replace_all(colnames(tmp), " ", "_"), "_avg_expr", sep = "")
  
  # pct
  for (i in levels(object)){
    tmp_cells <- Cells(subset(object, ident = i))
    
    tmp[,paste(str_replace_all(i, " ", "_"), "_pct", sep = "")] <- round(
      rowSums(x = object$ADT@data[rownames(deg_list), tmp_cells, drop = FALSE] > 0, na.rm = TRUE) /
        rowSums(!is.na(object$ADT@data[rownames(deg_list), tmp_cells, drop = FALSE])),
      digits = 3
    )
  }
  
  deg_list <- deg_list %>% 
    as_tibble(rownames = "Antibody")
  tmp <- tmp %>% 
    as_tibble(rownames = "Antibody")
  
  output_df <- deg_list %>% 
    left_join2(tmp)
  return(output_df)
}
