# Transcriptional Regulation of NP and NPR Signaling in Aging

This repository contains the code and results for studying the transcriptional regulation of neuropeptide (NP) and neuropeptide receptor (NPR) genes during aging in *Drosophila melanogaster*. The study explores how transcription factors (TFs) regulate these genes over time, with a focus on identifying TFs that control NP and NPR expression across different ages.

## Introduction

Aging involves significant gene expression changes controlled by transcription factors (TFs), especially those involved in the insulin/insulin-like growth factor 1 (IIS) pathway, the target of rapamycin (TOR) pathway, and the FOXO family of TFs. These pathways regulate aging processes such as antioxidant defense, DNA repair, and autophagy. We examine these processes in *Drosophila* using data from the Aging Fly Cell Atlas (AFCA).

### Key Findings
- NPR genes are regulated by more TFs than NP genes, suggesting a more complex regulatory network.
- At the mRNA level, NPR expression increases with age, while NP expression peaks at day 5 and declines thereafter.
- TF-NP pairs show a slight decrease in importance with age, while TF-NPR pairs have a more pronounced decline.
- The top-ranked TF-NPR pairs reveal stronger regulatory effects, indicating that NPR genes are influenced by key TFs with stronger effects.

The analysis uses scRNA-seq data from different ages (1, 3, 5, 20, 50, and 70 days) and employs SCENIC and other tools to predict TF regulons and importance scores.


## Code Overview

### 1. SCENIC Analysis
- **1_read_data.R**: Reads the raw data from the AFCA database.
- **2_for.scenic.data.R**: Processes the data for SCENIC analysis.
- **3_SCENIC.sh**: Shell script to run SCENIC analysis.

### 2. Network Analysis
- **importance_compare.R**: Compares the importance of TF-NP and TF-NPR regulatory pairs.

### 3. Visualization
- **heatmap.R**: Generates heatmaps to visualize the expression levels of TF-NP and TF-NPR pairs.
- **venn.R**: Creates Venn diagrams to show the overlap between TFs regulating NP and NPR genes.

## Data

The data acquisition process is described in the README.md file located in the data_aging folder.

## Results

The results of the analysis, including the importance scores and regulatory network diagrams, are stored in the `imgs` and `ave_importance` directories. These include PDFs of heatmaps, Venn diagrams, and network diagrams.

### `ave_importance` Folder

The `ave_importance` folder contains files related to the analysis of network importance for transcription factor (TF), neuropeptide (NP), and neuropeptide receptor (NPR) pairs. It is divided into two subfolders: `Top100_combined_TF-NP-NPR-pairs_network` and `Top50_seperated_TF-NP-NPR-pairs_network`. For a detailed description of the contents and usage of this folder, refer to the README file located within the `ave_importance` folder.

## Citation

Lu T-C, Brbić M, Park Y-J, Jackson T, Chen J, Kolluru SS, Qi Y, Katheder NS, Cai XT, Lee S, et al. 2023. Aging Fly Cell Atlas identifies exhaustive aging features at cellular resolution. Science. 380(6650):eadg0934. doi:10.1126/science.adg0934.

Aibar S, González-Blas CB, Moerman T, Huynh-Thu VA, Imrichova H, Hulselmans G, Rambow F, Marine J-C, Geurts P, Aerts J, et al. 2017. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 14(11):1083–1086. doi:10.1038/nmeth.4463.
