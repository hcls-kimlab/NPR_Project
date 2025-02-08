# Load required packages
library(igraph)
library(ggraph)
library(tidygraph)

#### Importance plot ####
# Read TF-gene data
# tf_gene_data <- read.csv("path_to_file.csv")
# Define file paths
files <- c("scenic_results/adj.afca_body_seurat_n_5.tsv",
           "scenic_results/adj.afca_head_seurat_n_5.tsv",
           "scenic_results/adj.afca_body_seurat_n_30.tsv",
           "scenic_results/adj.afca_head_seurat_n_30.tsv",
           "scenic_results/adj.afca_body_seurat_n_50.tsv",
           "scenic_results/adj.afca_head_seurat_n_50.tsv",
           "scenic_results/adj.afca_body_seurat_n_70.tsv",
           "scenic_results/adj.afca_head_seurat_n_70.tsv")

# Define NP and NPR lists
NPR_HR <- c("AkhR", "InR", "AstA-R1", "AstA-R2", "AstC-R1", "AstC-R2", "CapaR", "CCAP-R", 
            "CCHa1-R", "CCHa2-R", "CCKLR-17D1", "CCKLR-17D3", "CG4313", "CG12290", "CG13229", 
            "CG13575", "CG13995", "CG30340", "CG32547", "CG33639", "CNMaR", "CrzR", "ETHR", 
            "FMRFaR", "Lgr1", "Lgr3", "Lgr4", "Lkr", "moody", "MsR1", "MsR2", "NPFR", "PK1-R", 
            "PK2-R1", "PK2-R2", "Proc-R", "rk", "RYa-R", "SIFaR", "sNPF-R", "SPR", "TkR86C", 
            "TkR99D", "Tre1", "TrissinR")
NP <- c("Akh", "amn", "AstA", "AstC", "AstCC", "Burs", "Capa", "CCAP", "CCHa1", "CCHa2", 
        "CNMa", "Crz", "Dh31", "Dh44", "Dsk", "Eh", "ETH", "FMRFa", "Gpa2", "Gpb5", "Hug", 
        "Ilp1", "Ilp2", "Ilp3", "Ilp4", "Ilp5", "Ilp6", "Ilp7", "Ilp8", "ITP", "Lk", "Mip", 
        "Ms", "NPF", "Nplp1", "Nplp2", "Nplp3", "Nplp4", "Orcokinin", "Pburs", "Pdf", "Proc", 
        "Ptth", "RYa", "SIFa", "sNPF", "SP", "spab", "Tk")

# Read NPR list from file
NPR <- read.csv("id_validation_table_NPR-activity.txt", sep = "\t")
NPR <- NPR$current_symbol
NPR <- union(NPR, NPR_HR)

# Create a directory to save results
dir.create("ave_importance-organize", showWarnings = FALSE)

# Loop through each file
for (file in files) {
  
  # Read data
  tf_gene_data <- read.table(file, header = TRUE)
  
  # Filter NP and NPR target genes
  tf_gene_data_NP <- tf_gene_data[tf_gene_data$target %in% NP, ]
  tf_gene_data_NPR <- tf_gene_data[tf_gene_data$target %in% NPR, ]
  
  # Select top 50 NP and top 50 NPR
  tf_gene_data_NP <- tf_gene_data_NP[1:50, ]
  tf_gene_data_NPR <- tf_gene_data_NPR[1:50, ]
  
  # Merge selected NP and NPR data
  tf_gene_data_selected <- rbind(tf_gene_data_NP, tf_gene_data_NPR)
  
  # Add NP and NPR labels
  tf_gene_data_selected$gene_type <- ifelse(tf_gene_data_selected$target %in% NP, "NP", 
                                            ifelse(tf_gene_data_selected$target %in% NPR, "NPR", NA))
  
  # Calculate the average regulatory importance of each TF on NP and NPR
  importance_summary <- aggregate(importance ~ TF + gene_type, data = tf_gene_data_selected, FUN = mean)
  
  # Print result
  print(importance_summary)
  
  # Extract filename for saving
  file_name <- gsub("scenic_results/", "", file)
  file_name <- gsub(".tsv", ".importance_summary.csv", file_name)
  
  # Save result
  write.csv(importance_summary, paste0("ave_importance-organize/", file_name), row.names = FALSE)
}

### Set threshold to compare importance
for (file in files) {
  
  # Read data
  tf_gene_data <- read.table(file, header = TRUE)
  
  # Filter NP and NPR target genes
  tf_gene_data_NP <- tf_gene_data[tf_gene_data$target %in% NP, ]
  tf_gene_data_NPR <- tf_gene_data[tf_gene_data$target %in% NPR, ]
  
  # Merge selected NP and NPR data
  tf_gene_data_selected <- rbind(tf_gene_data_NP, tf_gene_data_NPR)
  tf_gene_data_selected <- tf_gene_data_selected[order(-as.numeric(tf_gene_data_selected$importance)),]
  
  # Select top 100 entries
  tf_gene_data_selected <- tf_gene_data_selected[1:100, ]
  
  # Add NP and NPR labels
  tf_gene_data_selected$gene_type <- ifelse(tf_gene_data_selected$target %in% NP, "NP", 
                                            ifelse(tf_gene_data_selected$target %in% NPR, "NPR", NA))
  
  # Calculate the average regulatory importance of each TF on NP and NPR
  importance_summary <- aggregate(importance ~ TF + gene_type, data = tf_gene_data_selected, FUN = mean)
  
  # Print result
  print(importance_summary)
  
  # Extract filename for saving
  file_name <- gsub("scenic_results/", "", file)
  file_name <- gsub(".tsv", ".importance_summary.csv", file_name)
  
  # Save result
  write.csv(importance_summary, paste0("ave_importance-organize/Top100_", file_name), row.names = FALSE)
}
