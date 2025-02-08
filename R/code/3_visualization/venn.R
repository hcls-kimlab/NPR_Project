## Visualization
# Clean the environment
rm(list=ls())

# Load necessary R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(loomR)
  library(SCopeLoomR)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(data.table)
  library(VennDiagram)
})

#### 1. Extract loom regulon information
path <- "D:/Windows/Share/FCAdata/"

# Read head data
loom_head <- Connect(filename = paste0(path, "s_fca_biohub_head_10x.loom"), mode = "r")
regulons_head <- regulonsToGeneLists(get_regulons(loom_head, column.attr.name="MotifRegulonGeneOccurrences"))
regulonAUC_head <- get_regulons_AUC(loom_head, column.attr.name='MotifRegulonsAUC')
close_loom(loom_head)

# Read body data
loom_body <- Connect(filename = paste0(path, "s_fca_biohub_body_10x.loom"), mode = "r")
regulons_body <- regulonsToGeneLists(get_regulons(loom_body, column.attr.name="MotifRegulonGeneOccurrences"))
regulonAUC_body <- get_regulons_AUC(loom_body, column.attr.name='MotifRegulonsAUC')
close_loom(loom_body)

#### 2. Process NP and NPR data
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

#### 3. Filter TF
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

tf_NP_head <- filter_tf_by_np_count(regulons_head, NP, threshold = 2)
tf_NPR_head <- filter_tf_by_np_count(regulons_head, NPR, threshold = 2)
tf_NP_body <- filter_tf_by_np_count(regulons_body, NP, threshold = 2)
tf_NPR_body <- filter_tf_by_np_count(regulons_body, NPR, threshold = 2)

tf_list <- list('tf_NP_head' = tf_NP_head, 'tf_NPR_head' = tf_NPR_head, 
                'tf_NP_body' = tf_NP_body, 'tf_NPR_body' = tf_NPR_body)

#### 4. Process TF results
tf_gene_list <- lapply(tf_list, function(tf_data) {
  sapply(names(tf_data), function(tf) str_sub(tf, end = -11))
})

data.frame(tf_gene_list$tf_NP_head, tf_gene_list$tf_NPR_head, tf_gene_list$tf_NP_body, tf_gene_list$tf_NPR_body)

#### 5. Venn diagram analysis
write.table(tf_gene_list$tf_NP_head, "Organized_results/Venn/tf_NP_head_gene_list.txt", row.names = FALSE, col.names = FALSE)
write.table(tf_gene_list$tf_NPR_head, "Organized_results/Venn/tf_NPR_head_gene_list.txt", row.names = FALSE, col.names = FALSE)
write.table(tf_gene_list$tf_NP_body, "Organized_results/Venn/tf_NP_body_gene_list.txt", row.names = FALSE, col.names = FALSE)
write.table(tf_gene_list$tf_NPR_body, "Organized_results/Venn/tf_NPR_body_gene_list.txt", row.names = FALSE, col.names = FALSE)

# Draw Venn diagram
venn1 <- venn.diagram(
  x = list(
    'NP Head' = tf_gene_list$tf_NP_head,
    'NPR Head' = tf_gene_list$tf_NPR_head,
    'NP Body' = tf_gene_list$tf_NP_body,
    'NPR Body' = tf_gene_list$tf_NPR_body
  ),
  scaled = TRUE, # Scale based on size
  alpha = 0.3, # Transparency
  lwd = 1, lty = 1, col = c('#85C4FC', '#FF9380','#0B8AF9','#FF2600'), # Line color
  label.col = 'black', # Label color
  cex = 1, # Label size
  fontface = "bold", # Bold font
  fill = c('#85C4FC', '#FF9380','#0B8AF9','#FF2600'), # Fill color
  category.names = c('NP Head', 'NPR Head', 'NP Body', 'NPR Body'), # Category names
  cat.dist = c(0.2, 0.2, 0.1, 0.1), # Category label distance
  cat.pos = c(-20, 20, -20, 20), # Category label angles
  cat.cex = 1, # Category label size
  cat.fontface = "bold", # Bold category labels
  cat.col = 'black', # Category label color
  resolution = 300, filename = NULL
)

# Save Venn diagram as PDF
pdf('Organized_results/Venn/VennDiagram.pdf', width = 3, height = 3)
grid.draw(venn1)
dev.off()
