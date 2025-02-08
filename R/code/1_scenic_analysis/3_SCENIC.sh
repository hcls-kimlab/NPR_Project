pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample1.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample1.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 3 \
--mask_dropouts

pyscenic aucell \
sample1.loom \
reg.csv \
--output sample1_SCENIC.loom \
--num_workers 3