# Data Sources

This project uses the following external datasets. Due to file size limitations, **please download them manually** and place them in the appropriate directories.

---

## 1. Cistarget Database

### Download Instructions
- **Source**: [Aertslab Cistarget Databases](https://resources.aertslab.org/cistarget/databases/)
- **Required Files**: 
  - `SCENIC motifs database` (`dm6_v10_clust.genes_vs_motifs.rankings.feather`)
- **Steps**:
  1. Visit the [download page](https://resources.aertslab.org/cistarget/databases/).
  2. Select the **"Fly (Drosophila melanogaster)"** databases.
  3. Download the files and place them in `data/cistarget/`.

## 2. Aging Fly Cell Atlas (AFCA)

### Download Instructions
- **Source**: [Hongjie Li Lab Web](https://hongjielilab.org/afca/)
- **Required Files**: 
  - `body data` (`adata_body_S_v1.0.h5ad`)
  - `head data` (`adata_head_S_v1.0.h5ad`)
- **Steps**:
  1. Visit the [download page](https://hongjielilab.shinyapps.io/AFCA/).
  2. Download the files and place them in `data/AFCA/raw`.

### Citation
```plaintext
Aibar, S., González-Blas, C., Moerman, T. et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods 14, 1083–1086 (2017).

Lu, T.-C., Brbić, M., Park, Y.-J., Jackson, T., Chen, J., Saroja, S., Qi, Y., Katheder, N.S., Cai, X.T., Lee, S., Chen, C., Auld, N., Liang, C.-Y., Ding, S.H., Welsch, D., Pisco, A.O., Jones, R.C., Leskovec, J., Lai, E.C., Luo, L., Jasper, H., Quake, S.R., Li, H., 2023. Aging Fly Cell Atlas Identifies Exhaustive Aging Features at Cellular Resolution.
