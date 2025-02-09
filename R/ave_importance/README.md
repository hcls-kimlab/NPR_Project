# `ave_importance` Folder Description

## Overview
The `ave_importance` folder contains files related to the analysis of network importance for transcription factor (TF), neuropeptide (NP), and neuropeptide receptor (NPR) pairs. It is divided into two subfolders: `Top100_combined_TF-NP-NPR-pairs_network` and `Top50_seperated_TF-NP-NPR-pairs_network`.

## Folder Structure

### 1. `Top100_combined_TF-NP-NPR-pairs_network`
This subfolder contains files related to the combined analysis of the Top100 TF-NP-NPR pairs.

- **`body_Top100.pdf`**: Network diagram showing the Top100 TF-NP-NPR pairs for the body region.
- **`head_Top100.pdf`**: Network diagram showing the Top100 TF-NP-NPR pairs for the head region.
- **`Top100.prism`**: Data file for analysis using Prism software.
- **`Top100_adj.afca_body_seurat_n_*.importance_summary.csv`**: Importance summary data for the body region under different parameters (n=5, 30, 50, 70).
- **`Top100_adj.afca_head_seurat_n_*.importance_summary.csv`**: Importance summary data for the head region under different parameters (n=5, 30, 50, 70).

### 2. `Top50_seperated_TF-NP-NPR-pairs_network`
This subfolder contains files related to the separated analysis of the Top50 TF-NP-NPR pairs.

- **`adj.afca_body_seurat_n_*.importance_summary.csv`**: Importance summary data for the body region under different parameters (n=5, 30, 50, 70).
- **`adj.afca_head_seurat_n_*.importance_summary.csv`**: Importance summary data for the head region under different parameters (n=5, 30, 50, 70).
- **`ave_importance.prism`**: Data file for analysis using Prism software.
- **`body.pdf`**: Network diagram showing the Top50 TF-NP-NPR pairs for the body region.
- **`head.pdf`**: Network diagram showing the Top50 TF-NP-NPR pairs for the head region.
- **`Top50_seperated_TF-NP-NPR-pairs_network.pdf`**: Separated network diagram of the Top50 TF-NP-NPR pairs.

## Usage Instructions
- Use the **`Top100.prism`** and **`ave_importance.prism`** files for further statistical analysis in Prism.
- The **`body_Top100.pdf`**, **`head_Top100.pdf`**, **`body.pdf`**, and **`head.pdf`** files can be used to visualize the network structures.
- The **`importance_summary.csv`** files provide detailed importance scores for further analysis.

## Notes
- Ensure that the latest version of Prism software is installed when working with `.prism` files.
- All CSV files can be opened and viewed using Excel or any text editor.

For any questions, please contact the project lead.
