import loompy
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import scipy

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns

filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/graph/gpcr"
filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/h5ad"
filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_body_10x.loom"
filename="s_fca_biohub_head_10x.loom"
#filename="Gr5a_leg_10x.loom"
filename="Abdominal_ganglion_cluster_circle.loom"
filename="s_fca_biohub_haltere_10x.loom"
#filename ="r_fca_biohub_testis_10x.loom"

partid = 8
#filename="s_fca_biohub_head_10x.h5ad"

if partid==3:
    clockpath="D:/project/python/bioinfo/data/graph/clock/body/"
elif partid==8:
    clockpath = "D:/project/python/bioinfo/data/graph/clock/head/"


#filename="s_fca_biohub_head_10x.h5ad"
#filename="s_fca_biohub_head_10x.loom"
#filename="s_fca_biohub_head_10x_.loom"
#filename="fca_biohub_fat_body_ss2.loom"
#filename="s_fca_biohub_head_10x.h5ad"
#filename="s_fca_biohub_heart_10x.loom"

filedir= os.path.join(filepath,filename)

print(filedir)

#with loompy.connect(filedir) as ds:
#    print(ds)

#adata = anndata.AnnData

#print("===========adata=================")
#print(adata)
#adata.X =

ds = loompy.connect(filedir)
print("===========ds=================")
#12158 6527
print(ds.shape)
print(ds)


nparr = np.array(ds[:, :], dtype='int')
print("=============ds nparr======================")
print(nparr)
anndataX= scipy.sparse.csc_matrix(nparr)

#scipy.io.mmwrite(os.path.join(pheromone,mfilename), narry)

# load sparse matrix:
#X = io.mmread("counts.mtx")
# create anndata object
#adata = anndata.AnnData(
#    X=X.transpose().tocsr()
#)


print("=============h5ad anndataadata.shape======================")

print(anndataX.shape)

#0,0 로 나온다 하지만 mtx로 변경할때 1,1로 넣어줘야 한다..


#totals = ds.map([np.sum], axis=1)[0]  # Calculate the total molecule count for each cell
#print("total==>",totals)
#cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
#print("cells==>",len(cells))


print("=====col_attrs==============")
print(ds.col_attrs)
print("=====row_attrs==============")
print(ds.row_attrs)
#totals = ds.map([np.sum], axis=1)[0]
#cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
print("========tissue===========")
print(ds.col_attrs['tissue'])

print("=======sex==============")
sex= np.array(ds.col_attrs['sex'])
print(sex.shape)
print(sex)

print("=======sex==============")
cell= np.array(ds.col_attrs['CellID'])
print(cell.shape)
print(cell)

print("=======age==============")
age= np.array(ds.col_attrs['age'])
print(age)
print(len(age))

#print(ds.row_attrs)
print("=======gene==============")
gene= np.array(ds.row_attrs['Gene'])

print(gene)
print(len(gene))
print(gene[0].decode('utf8'))



print("===========barcode==============")
print(ds.col_attrs['Barcode'])
print(len(ds.col_attrs['Barcode']))

print("===========annotation__ontology_id==============")
print(ds.col_attrs['annotation__ontology_id'])
print(len(ds.col_attrs['annotation__ontology_id']))

print("===========cell id==============")
print(ds.col_attrs['CellID'])
print(len(ds.col_attrs['CellID']))

print("===========annontation===========")
print(ds.col_attrs['annotation'])
print(len(ds.col_attrs['annotation']))  #유전자별로..

print("===========MotifRegulons================")
print(ds.row_attrs['MotifRegulons'])
print(len(ds.row_attrs['MotifRegulons']))  #유전자별로..

print(ds.row_attrs['MotifRegulons'])
print(len(ds.row_attrs['MotifRegulons']))  #유전자별로..

print("===========MotifRegulons= Abd-B_(+)-motif MotifRegulons===============")
print(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'])  #유전자별로..
print(ds.row_attrs['MotifRegulons'].dtype)  #314개가 존재한다.

print("===========MotifRegulons dtype==============")
#print(ds.row_attrs['MotifRegulons'].dtype)
#tf = np.array(ds.row_attrs['MotifRegulons'].dtype)

tf = [x for x,y in ds.row_attrs['MotifRegulons'].dtype.fields.items()]
print(len(tf))
print(tf[0])  #motif full name
print(tf[0].index('_('))
print(tf[0][:tf[0].index('_(')])  #motif gene name
#'MotifRegulonGeneWeights
print(ds.row_attrs['MotifRegulons'][tf[0]]) #matrix를 가져온다..

#순서
print(ds.row_attrs['MotifRegulons'][tf[0]][20]) #matrix를 가져온다..
print(ds.row_attrs['MotifRegulonGeneWeights'][tf[0]][20]) #matrix를 가져온다..
print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'][12094])

motifregeulons = np.transpose(np.nonzero(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif']))
print(motifregeulons)



motifregeulons1=np.transpose(np.nonzero(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif']))
motifregeulonsw= np.transpose(np.nonzero(ds.row_attrs['MotifRegulonGeneWeights']['Abd-B_(+)-motif']))
#print(motifregeulons1)
print(len(motifregeulons1))
print(len(motifregeulonsw))
#print(motifregeulons[20])
MotifRegulonGeneOccurrences = np.transpose(np.nonzero(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif']))
print(len(MotifRegulonGeneOccurrences))
print(ds.row_attrs['MotifRegulonGeneOccurrences'].shape)  #9076
print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'])

print(motifregeulons1[0][0])
print(motifregeulonsw[0][0])
print(MotifRegulonGeneOccurrences[0][0])


print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'][12094])
print(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'][12094]) #matrix를 가져온다..
print(ds.row_attrs['MotifRegulonGeneWeights']['Abd-B_(+)-motif'][12094]) #matrix를 가져온다..

print(ds.col_attrs['MotifRegulonsAUC'].dtype)
print(ds.col_attrs['MotifRegulonsAUC'].shape)

tfauc = [x for x,y in ds.col_attrs['MotifRegulonsAUC'].dtype.fields.items()]
print(len(tfauc)) #314
MotifRegulontfauc = np.transpose(np.nonzero(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif']))
print(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif'][0])
#6517개가 반복된다..
print(len(MotifRegulontfauc))
print(len(ds.col_attrs['MotifRegulonsAUC']['Abd-B_(+)-motif']))
print(len(ds.col_attrs['MotifRegulonsAUC']))
#6527


'''
for regrow in motifregeulons:
    print(regrow[0])

    for generow in motifregeulons:
        print(generow[0])  # gene matrix
        print(ds.row_attrs['MotifRegulonGeneOccurrences']['Abd-B_(+)-motif'][generow[0]])
        print(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'][generow[0]])  # matrix를 가져온다..
        print(ds.row_attrs['MotifRegulonGeneWeights']['Abd-B_(+)-motif'][generow[0]])  # matrix를 가져온다..

'''

#MotifRegulonGeneOccurrences
#MotifRegulonGeneWeights
#MotifRegulons
'''
for i in np.transpose(np.nonzero(ds.row_attrs['MotifRegulons']['Abd-B_(+)-motif'])):
    print(i[0])
    print(gene[i[0]].decode('utf8'))  #genename->
'''

#314개..

#adata =anndata.AnnData(anndataX.T,obs=ds.col_attrs , var = ds.row_attr)

#obs_names = "cell_names",
#var_names = "gene_names",


print("================adata var===============")
#print(adata.var)

#Abd-B_(+)-motif
'''
get_regulons(
  loom,
  column.attr.name = "MotifRegulons",
  tf.as.name = TRUE,
  tf.sep = "_"
)
# Regulons:
regulons_incidMat <- get_regulons(loom) # as incid matrix
regulons_motif <- regulonsToGeneLists(regulons_incidMat) # convert to list
regulons_ChIP <- regulonsToGeneLists(get_regulons(loom, attrName = "TrackRegulons"))
# Regulon AUC and thresholds
regulonsAUC <- get_regulonsAuc(loom)
regulonsAucThresholds <- get_regulonThresholds(loom)
# Embeddings (tsne/umap)
embeddings <- get_embeddings(loom)
'''
print("================adata obs===============")
#print(adata.obs)

ds.close()