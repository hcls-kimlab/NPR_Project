TFs <- read.csv("cisTarget_databases/allTFs_dmel.txt", header = F)$V1

## Get the NP and NPR list
NPR_HR <- c("AkhR","InR", "AstA-R1", "AstA-R2", "AstC-R1", "AstC-R2", "CapaR", "CCAP-R", "CCHa1-R", "CCHa2-R", "CCKLR-17D1", "CCKLR-17D3", "CG4313", "CG12290", "CG13229", "CG13575", "CG13995", "CG30340", "CG32547", "CG33639", "CNMaR", "CrzR", "ETHR", "FMRFaR", "Lgr1", "Lgr3", "Lgr4", "Lkr", "moody", "MsR1", "MsR2", "NPFR", "PK1-R", "PK2-R1", "PK2-R2", "Proc-R", "rk", "RYa-R", "SIFaR", "sNPF-R", "SPR", "TkR86C", "TkR99D", "Tre1", "TrissinR")
NP <- c("Akh", "amn", "AstA", "AstC", "AstCC", "Burs", "Capa", "CCAP", "CCHa1", "CCHa2", "CNMa", "Crz", "Dh31", "Dh44", "Dsk", "Eh", "ETH", "FMRFa", "Gpa2", "Gpb5", "Hug", "Ilp1", "Ilp2", "Ilp3", "Ilp4", "Ilp5", "Ilp6", "Ilp7", "Ilp8", "ITP", "Lk", "Mip", "Ms", "NPF", "Nplp1", "Nplp2", "Nplp3", "Nplp4", "Orcokinin", "Pburs", "Pdf", "Proc", "Ptth", "RYa", "SIFa", "sNPF", "SP", "spab", "Tk")

# Load NPR data from a file
NPR <- read.csv("id_validation_table_NPR-activity.txt", sep = "\t")
NPR <- NPR$current_symbol
NPR <- union(NPR, NPR_HR)

# Combine TFs, NP, and NPR gene lists
genes_TF_NP_NPR <- union(TFs, NP)
genes_TF_NP_NPR <- union(genes_TF_NP_NPR, NPR)

# Check the first few entries
genes_TF_NP_NPR

# Subset Seurat objects to include only relevant genes
afca_head_seurat_n <- subset(afca_head_seurat, features = genes_TF_NP_NPR)
afca_body_seurat_n <- subset(afca_body_seurat, features = genes_TF_NP_NPR)

# Identify cells for different ages in the head Seurat object
Cells.head.5 <- WhichCells(afca_head_seurat_n, expression = age == '5')
Cells.head.30 <- WhichCells(afca_head_seurat_n, expression = age == '30')
Cells.head.50 <- WhichCells(afca_head_seurat_n, expression = age == '50')
Cells.head.70 <- WhichCells(afca_head_seurat_n, expression = age == '70')

# Identify cells for different ages in the body Seurat object
Cells.body.5 <- WhichCells(afca_body_seurat_n, expression = age == '5')
Cells.body.30 <- WhichCells(afca_body_seurat_n, expression = age == '30')
Cells.body.50 <- WhichCells(afca_body_seurat_n, expression = age == '50')
Cells.body.70 <- WhichCells(afca_body_seurat_n, expression = age == '70')

# Write count data for different cell groups in the head Seurat object
write.csv(t(as.matrix(afca_head_seurat_n[,Cells.head.5]@assays$RNA@counts)), file = "afca_head_seurat_n_5.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[,Cells.head.30]@assays$RNA@counts)), file = "afca_head_seurat_n_30.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[,Cells.head.50]@assays$RNA@counts)), file = "afca_head_seurat_n_50.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[,Cells.head.70]@assays$RNA@counts)), file = "afca_head_seurat_n_70.for.scenic.data.csv")

# Write count data for different cell groups in the body Seurat object
write.csv(t(as.matrix(afca_body_seurat_n[,Cells.body.5]@assays$RNA@counts)), file = "afca_body_seurat_n_5.for.scenic.data.csv") 
write.csv(t(as.matrix(afca_body_seurat_n[,Cells.body.30]@assays$RNA@counts)), file = "afca_body_seurat_n_30.for.scenic.data.csv")
write.csv(t(as.matrix(afca_body_seurat_n[,Cells.body.50]@assays$RNA@counts)), file = "afca_body_seurat_n_50.for.scenic.data.csv")
write.csv(t(as.matrix(afca_body_seurat_n[,Cells.body.70]@assays$RNA@counts)), file = "afca_body_seurat_n_70.for.scenic.data.csv")

# Optionally, uncomment to write all body Seurat object count data
# write.csv(t(as.matrix(afca_body_seurat_n@assays$RNA@counts)), file = "afca_body_seurat_n.for.scenic.data.csv")
