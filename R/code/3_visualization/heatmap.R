# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pheatmap)
library(stringr)

# Function to filter TFs based on the count of NP genes
filter_tf_by_np_count <- function(regulons, NP_list, threshold = 5) {
  result <- list()
  for (tf in names(regulons)) {
    np_count <- sum(regulons[[tf]] %in% NP_list)
    if (np_count >= threshold) {
      result[[tf]] <- np_count
    }
  }
  return(result)
}

# Get NP and NPR lists
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
NPR <- read.csv("id_validation_table_NPR-activity.txt", sep = "\t")$current_symbol
NPR <- union(NPR, NPR_HR)

# 1. Extract loom regulon information
path = "D:/Windows/Share/FCAdata/"

# Load loom data for 'head' and 'body'
loom_head <- Connect(filename = paste0(path, "s_fca_biohub_head_10x.loom"), mode = "r")
loom_body <- Connect(filename = paste0(path, "s_fca_biohub_body_10x.loom"), mode = "r")

# Extract regulon matrices and other necessary information
regulons_incidMat_head <- get_regulons(loom_head, column.attr.name = "MotifRegulonGeneOccurrences")
regulons_head <- regulonsToGeneLists(regulons_incidMat_head)
regulonAUC_head <- get_regulons_AUC(loom_head, column.attr.name = 'MotifRegulonsAUC')
regulonAucThresholds_head <- get_regulon_thresholds(loom_head)

regulons_incidMat_body <- get_regulons(loom_body, column.attr.name = "MotifRegulonGeneOccurrences")
regulons_body <- regulonsToGeneLists(regulons_incidMat_body)
regulonAUC_body <- get_regulons_AUC(loom_body, column.attr.name = 'MotifRegulonsAUC')
regulonAucThresholds_body <- get_regulon_thresholds(loom_body)

# Close loom files
close_loom(loom_head)
close_loom(loom_body)

# Filter TFs based on NP and NPR count for both head and body
tf_NP_head <- filter_tf_by_np_count(regulons_head, NP, threshold = 2)
tf_NPR_head <- filter_tf_by_np_count(regulons_head, NPR, threshold = 2)
tf_NP_body <- filter_tf_by_np_count(regulons_body, NP, threshold = 2)
tf_NPR_body <- filter_tf_by_np_count(regulons_body, NPR, threshold = 2)

# Combine NP and NPR TFs into one list
tf_select <- filter_tf_by_np_count(regulons_body, c(NPR, NP), threshold = 2)

# Process TF lists and generate gene lists for further analysis
tf_list <- list('tf_NP_head' = tf_NP_head, 'tf_NPR_head' = tf_NPR_head, 
                'tf_NP_body' = tf_NP_body, 'tf_NPR_body' = tf_NPR_body)
tf_gene_list <- list()

for (j in c('tf_NP_head', 'tf_NPR_head', 'tf_NP_body', 'tf_NPR_body')) {
  x <- c()
  y <- c()
  for (i in names(tf_list[[j]])) {
    x <- c(x, paste0(str_sub(i, end = -10), str_extract(j, "[NPR]+"), "_", as.character(tf_list[[j]][i]), "_regulon_ave"))
    y <- c(y, str_sub(i, end = -11))
  }
  tf_list[[j]] <- x
  tf_gene_list[[j]] <- y
}

# Output TF lists and gene lists
tf_list
tf_gene_list

# Body expression data
exp_NP_body <- AverageExpression(afca_body_seurat, group.by = c('age'), features = NP)$RNA
exp_NP_body %>% colMeans(na.rm = TRUE) %>% as.data.frame() %>% t() -> exp_NP_body
rownames(exp_NP_body) <- "NP_ave"

exp_NPR_body <- AverageExpression(afca_body_seurat, group.by = c('age'), features = NPR)$RNA
exp_NPR_body %>% colMeans(na.rm = TRUE) %>% as.data.frame() %>% t() -> exp_NPR_body
rownames(exp_NPR_body) <- "NPR_ave"

# Get expression data for TF genes
exp_tf_NP_body <- AverageExpression(afca_body_seurat, group.by = c('age'), features = tf_gene_list$tf_NP_body)$RNA
exp_tf_NPR_body <- AverageExpression(afca_body_seurat, group.by = c('age'), features = tf_gene_list$tf_NPR_body)$RNA

# Combine NP and NPR expression data for body
exp_body_NP <- rbind(exp_NP_body, exp_tf_NP_body)
exp_body_NP <- as.data.frame(exp_body_NP)
exp_body_NP[, "tf_type"] <- "NP"

exp_body_NPR <- rbind(exp_NPR_body, exp_tf_NPR_body)
exp_body_NPR <- as.data.frame(exp_body_NPR)
exp_body_NPR[, "tf_type"] <- "NPR"

# Remove overlapping rows and add a 'both' label
tf_both <- intersect(rownames(exp_body_NP), rownames(exp_body_NPR))
df_body <- rbind(exp_body_NP[!rownames(exp_body_NP) %in% rownames(exp_body_NPR), ], exp_body_NPR)
df_body[tf_both, ]$tf_type <- 'both'

# Head expression data
exp_NP_head <- AverageExpression(afca_head_seurat, group.by = c('age'), features = NP)$RNA
exp_NP_head %>% colMeans(na.rm = TRUE) %>% as.data.frame() %>% t() -> exp_NP_head
rownames(exp_NP_head) <- "NP_ave"

exp_NPR_head <- AverageExpression(afca_head_seurat, group.by = c('age'), features = NPR)$RNA
exp_NPR_head %>% colMeans(na.rm = TRUE) %>% as.data.frame() %>% t() -> exp_NPR_head
rownames(exp_NPR_head) <- "NPR_ave"

# Get expression data for TF genes
exp_tf_NP_head <- AverageExpression(afca_head_seurat, group.by = c('age'), features = tf_gene_list$tf_NP_head)$RNA
exp_tf_NPR_head <- AverageExpression(afca_head_seurat, group.by = c('age'), features = tf_gene_list$tf_NPR_head)$RNA

# Combine NP and NPR expression data for head
exp_head_NP <- rbind(exp_NP_head, exp_tf_NP_head)
exp_head_NP <- as.data.frame(exp_head_NP)
exp_head_NP[, "tf_type"] <- "NP"

exp_head_NPR <- rbind(exp_NPR_head, exp_tf_NPR_head)
exp_head_NPR <- as.data.frame(exp_head_NPR)
exp_head_NPR[, "tf_type"] <- "NPR"

# Remove overlapping rows and add a 'both' label
tf_both <- intersect(rownames(exp_head_NP), rownames(exp_head_NPR))
df_head <- rbind(exp_head_NP[!rownames(exp_head_NP) %in% rownames(exp_head_NPR), ], exp_head_NPR)
df_head[tf_both, ]$tf_type <- 'both'

# Heatmap plotting (body data)
df_body_sorted <- df_body[order(df_body$tf_type), ]
p <- pheatmap(as.matrix(df_body_sorted[, -ncol(df_body_sorted)]), 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,annotation_colors = ann_colors,
         scale = "row", color = colorRampPalette(brewer.pal(9, "GnBu"))(100))
p$tree_col[]$order # Order
# Heatmap plotting (head data)
df_head_sorted <- df_head[order(df_head$tf_type), ]
p <- pheatmap(as.matrix(df_head_sorted[, -ncol(df_head_sorted)]), 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,annotation_colors = ann_colors,
         scale = "row", color = colorRampPalette(brewer.pal(9, "GnBu"))(100))
p$tree_col[]$order # Order
