#np , npr , tf drawgraph each tissue in one figure

import loompy
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import desc
import scipy
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns

import matplotlib.font_manager as fm
import networkx as nx

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib.pyplot import rc_context

#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns


#tf and np network
#size up node display different

filepath="D:/project/python/bioinfo/data/loom"
#filepath="D:/project/python/bioinfo/data/h5ad"
filepath="D:/project/python/bioinfo/data/graph"
hg19="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/hg19/"
chgnew="D:/project/python/bioinfo/data/graph/filtered_gene_bc_matrices/chgnew"
pheromone="D:/project/python/bioinfo/data/graph/creatematrix/pheromone/"
creatematrixpath="D:/project/python/bioinfo/data/graph/creatematrix"
chemopath="D:/project/python/bioinfo/data/graph/chemo"
tfpath="D:/project/python/bioinfo/data/graph/tf/except"

filename="r_fca_biohub_testis_10x.loom"
filename="s_fca_biohub_antenna_10x.loom"
filename="s_fca_biohub_head_10x.loom"
filename="Chemosensory.txt"
filename="np_10_regulon_tissue_ann_tfgene_reggene_list.txt"
filename="figure7_np+npr_np_npr_10_no_h_b_all_network_npr.txt"   #np tissue 10
filename="figure7_np_npr_not_exists_whole.txt"
savepath="D:/project/python/bioinfo/figures/tf/network/except"

gfilename="genes.tsv"
bfilename="barcodes.tsv"
mfilename="matrix.mtx"
annofilename="pheromone_ann.tsv"
#filename="s_fca_biohub_head_10x.h5ad"
#filename="s_fca_biohub_head_10x.loom"
#filename="s_fca_biohub_head_10x_.loom"
#filename="fca_biohub_fat_body_ss2.loom"
#filename="s_fca_biohub_head_10x.h5ad"
#filename="10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5"
results_dir = 'D:/project/python/bioinfo/data/h5ad'  # the file that will store the analysis results
results_file="graph_5gene.h5ad"

#For whole tissue

filedir= os.path.join(tfpath,filename)

print(filedir)

#pbarcodes = pd.DataFrame(os.path.join(chemopath,filename),index=0)
dtfm = pd.read_csv(os.path.join(tfpath,filename), delimiter='\t')
#tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight

print(dtfm)
print(dtfm['tname'])
print(dtfm)

dtnp = dtfm[dtfm['nprtype'] == 'np']
dtnpr = dtfm[dtfm['nprtype'] == 'npr']

nprtdata = dtfm

print(dtnp)  #np
print(dtnpr) #npr




'''
tname_annotation = pd.DataFrame(
nprtdata.groupby(['tname', 'annotation','tfgenename','rgenename'])["rgenename"].count()).rename(columns={"rgenename": "Count"})
print("==========tname_annotation=======")
print(tname_annotation)
'''

ltname=['Testis','Antenna','Body wall','Fat body','Gut','Haltere','Heart','Leg','Male reproductive glands','Malpighian tubule',
        'Oenocyte','Ovary','Proboscis','Trachea','Wing','Head','Body']


def drawGraph(ltname):

    for j in range(17):

        print("===i===",j)
        stname = ltname[j]

        dtnp = dtfm[(dtfm['nprtype'] == 'np') & (dtfm['tname'] == stname)]
        dtnpr = dtfm[(dtfm['nprtype'] == 'npr') & (dtfm['tname'] == stname)]

        nprtdata = dtfm

        print(dtnp)  # np
        print(dtnpr)  # npr


        if(len(dtnp) >0):

            nptfgenenamedict = {}
            npuniqtfgenename = np.unique(dtnp.tf)
            print("======npuniqtfgenename[0]=========", npuniqtfgenename[0])

            for i in range(0, len(np.unique(dtnp.tf))):
                nptfgenenamedict[npuniqtfgenename[i]] = i

            npgenenamedict = {}
            npuniqrgenename = np.unique(dtnp.gene)
            print("======agenename[0]=========", npuniqrgenename[0])
            for i in range(0, len(np.unique(dtnp.gene))):
                npgenenamedict[npuniqrgenename[i]] = i

        if (len(dtnpr) > 0):
            nprtfgenenamedict = {}
            npruniqtfgenename = np.unique(dtnpr.tf)
            print("======agenename[0]=========", npruniqtfgenename[0])
            for i in range(0, len(np.unique(dtnpr.tf))):
                nprtfgenenamedict[npruniqtfgenename[i]] = i

            nprgenenamedict = {}
            npruniqrgenename = np.unique(dtnpr.gene)
            print("======agenename[0]=========", npruniqrgenename[0])
            for i in range(0, len(np.unique(dtnpr.gene))):
                nprgenenamedict[npruniqrgenename[i]] = i

            print("=======nprgenenamedict=====")
            print(nprgenenamedict)
            print("=======npruniqrgenename=====")
            print(npruniqrgenename)

        # tname	motifname	tfgenename	tfcount	npgenename	npmotifauc	npgeneexpress	npmotifregoccur	npmotifregweight


        # npgene sum(count)
        if (len(dtnp) > 0):
            npgenecountsum = dtnp.groupby('gene')['gcount'].sum()
            npmotifnamesum = dtnp.groupby('tf')['tfcount'].sum()

        if (len(dtnpr) > 0):
            nprgenecountsum = dtnpr.groupby('gene')['gcount'].sum()
            nprmotifnamesum = dtnpr.groupby('tf')['tfcount'].sum()


        G = nx.Graph(day="Stackoverflow")

        for index, row in dtnp.iterrows():
            # G.add_node(row['motifname'], group=row['npgenename'], nodesize=row['tfcount'])
            G.add_node(row['tf'], nodesize=1)

        for index, row in dtnp.iterrows():
            # G.add_weighted_edges_from([(row['motifname'], row['npgenename'], row['npmotifregweight']*100)])
            G.add_node(row['gene'], nodesize=1)

        for index, row in dtnp.iterrows():
            G.add_edges_from([(row['tf'], row['gene'])])

        for index, row in dtnpr.iterrows():
            # G.add_node(row['motifname'], group=row['npgenename'], nodesize=row['tfcount'])
            G.add_node(row['tf'], nodesize=1)

        for index, row in dtnpr.iterrows():
            # G.add_weighted_edges_from([(row['motifname'], row['npgenename'], row['npmotifregweight']*100)])
            G.add_node(row['gene'], nodesize=1)

        for index, row in dtnpr.iterrows():
            G.add_edges_from([(row['tf'], row['gene'])])

        color_map = {1: '#f09494', 2: '#eebcbc', 3: '#72bbd0', 4: '#91f0a1', 5: '#629fff', 6: '#bcc2f2',
                     7: '#eebcbc', 8: '#f1f0c0', 9: '#d2ffe7', 10: '#caf3a6', 11: '#ffdf55', 12: '#ef77aa',
                     13: '#d6dcff', 14: '#d2f5f0'}
        # tf	tfcount	gene	gcount	nprtype
        # dtnp , dtnpr
        # nptfgenenamedict, npuniqtfgenename ,npgenenamedict ,npuniqrgenename
        # nprtfgenenamedict,npruniqtfgenename,nprgenenamedict,npruniqrgenename

        colors = []
        nptfgenecount = dtnp.tf.value_counts()
        npgenecount = dtnp.gene.value_counts()

        nprtfgenecount = dtnpr.tf.value_counts()
        nprgenecount = dtnpr.gene.value_counts()

        print("===Bodyrgenecountsum===")

        node_sizes = []
        for node in G:
            # 6F495C  tf npr
            # B791A4  tf np
            # 72bbd0 npr
            # caf3a6 np
            nptfgenecount = dtnp.tf.value_counts()
            npgenecount = dtnp.gene.value_counts()

            nprtfgenecount = dtnpr.tf.value_counts()
            nprgenecount = dtnpr.gene.value_counts()

            if (nptfgenecount.get(node) != None):
                colors.append('#000000')  # 검은색
            elif (npgenecount.get(node) != None):
                colors.append('#8C8C8C')  # 회색
            elif (nprtfgenecount.get(node) != None):
                colors.append('#FF0000')  # 빨간색
            elif (nprgenecount.get(node) != None):
                colors.append('#FFBB00')  # 오렌지

        # colors = [color_map[G.node[node]['rgenename']] for node in G]
        # sizes = [G.node[node]['nodesize']*10 for node in G]
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['font.family'] = 'Helvetica'  # 'Helvetica'
        plt.rcParams['font.size'] = 5

        plt.figure(figsize=(13, 13))
        plt.title(stname, fontsize=20)

        options = {
            'edge_color': '#FFDEA2',
            'node_size': 3,
            'width': 1,
            'with_labels': True,
            'font_weight': 'regular',
            'font_family': 'Helvetica',
            'font_size': 5,
        }

        f = plt.figure(1)
        # ax = f.add_subplot(1,1,1)
        ax = plt.gca()

        # print('====rgenecountsum.index(0)====')
        # print(rgenecountsum.index[0])
        # print(rgenecountsum.get(0))

        # print(motifnamesum.index[0])
        # print(motifnamesum.get(0))
        # ('He is %s %s and is %i years old.') % (first_name, last_name, age)


        legend_elements = [Line2D([0], [0], color='#FFDEA2', label='Edge'),
                           Line2D([0], [0], marker='o', markersize=10, color='#000000',
                                  label='TF-NP'),
                           Line2D([0], [0], marker='o', markersize=10, color='#FF0000',
                                  label='TF-NPR'),
                           Line2D([0], [0], marker='o', markersize=10, color='#8C8C8C',
                                  label='NP'),
                           Line2D([0], [0], marker='o', markersize=10, color='#FFBB00', label='NPR')]

        # nx.draw(G, node_color=colors, node_size=sizes,pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
        # nx.draw(G, node_color=colors,  pos=nx.spring_layout(G, k=1.0, iterations=6), **options)
        nx.draw_networkx(G, node_color=colors, pos=nx.spring_layout(G, k=0.25, iterations=50), **options, ax=ax)

        f.set_facecolor('w')
        plt.legend(handles=legend_elements, fontsize=12, loc='upper left')
        f.tight_layout()

        savepath = "D:/project/python/bioinfo/figures/tf/network/except"
        savefile = "figure7_tf_np_npr_network_" + ltname[j] + "_except.pdf"
        plt.savefig(os.path.join(savepath, savefile))

        savepath = "D:/project/python/bioinfo/figures/tf/network/except"
        savefile = "figure7_tf_np_npr_network_" + ltname[j] + "_except.png"
        plt.savefig(os.path.join(savepath, savefile))

        plt.show()
        print(ltname[j])


#each tissue shows tf, np
drawGraph(ltname)




