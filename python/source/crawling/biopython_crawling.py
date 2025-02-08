import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import numpy as np
import logging
import os
import sys
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
import requests
import json


print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))  #D:\data\slb\project\python
basedir = os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))) #D:\data\slb\project\python
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(__file__))) #load_info

import sysinfo
import platform
import pymysql
'''
Aasdh	FBgn0027780
Ac3	FBgn0023416
Ac76E	FBgn0004852
AcCoAS	FBgn0012034
Ace	FBgn0000024
Acn	FBgn0263198
Act5C	FBgn0000042
Adar	FBgn0026086
Adk1	FBgn0037995
Ae2	FBgn0036043
'''
url = "http://flybase.org/api/sequence/id/FBgn0025595"
resp=requests.get(url)
jstring = json.dumps(resp.text)

lgene=[
'AkhR',
'AstA-R1',
'AstA-R2',

'AstC-R1',
'AstC-R2',
'rk',
'rk',
'CapaR',
'PK1-R',
'CCHa1-R',
'CCHa2-R',
'CNMaR',
'CrzR',
'CCAP-R',
'Dh31-R',

'Dh44-R1',
'Dh44-R2',
'ETHR',

'FMRFaR',
'Lgr1',
'Lgr1',
'PK2-R2',
'PK2-R1',
'InR',
'InR',
'InR',
'InR',
'InR',
'InR',
'Lgr4',
'Lgr3',
'Lkr',
'PK1-R',
'MsR1',
'MsR2',
'TkR86C',
'NPFR',
'Gyc76C',
'Pdfr',
'Proc-R',

'RYa-R',
'SPR',
'sNPF-R',
'SIFaR',
'CCKLR-17D1',
'CCKLR-17D3',
'TkR99D',
'TrissinR'



]

lsys=[
'FBgn0025595',
'FBgn0266429',
'FBgn0039595',

'FBgn0036790',
'FBgn0036789',
'FBgn0003255',
'FBgn0003255',
'FBgn0037100',
'FBgn0038201',
'FBgn0050106',
'FBgn0033058',
'FBgn0053696',
'FBgn0036278',
'FBgn0039396',
'FBgn0052843',

'FBgn0033932',
'FBgn0033744',
'FBgn0038874',

'FBgn0035385',
'FBgn0016650',
'FBgn0016650',
'FBgn0038139',
'FBgn0038140',
'FBgn0283499',
'FBgn0283499',
'FBgn0283499',
'FBgn0283499',
'FBgn0283499',
'FBgn0283499',
'FBgn0085440',
'FBgn0039354',
'FBgn0035610',
'FBgn0038201',
'FBgn0035331',
'FBgn0264002',
'FBgn0004841',
'FBgn0037408',
'FBgn0266136',
'FBgn0260753',
'FBgn0029723',

'FBgn0004842',
'FBgn0029768',
'FBgn0036934',
'FBgn0038880',
'FBgn0259231',
'FBgn0030954',
'FBgn0004622',
'FBgn0085410'


]
print("=========jstring========")
#print(jstring)

arr = jstring.split('length=')
print(arr[1])
lstring = arr[1]
arrl = lstring.split(';')
print("===========length==========")
print(arrl[0])

print("=======lsys========")
for i in range(len(lsys)):
    #print(lsys[0])
    url = "http://flybase.org/api/sequence/id/"+lsys[i]
    resp = requests.get(url)
    jstring = json.dumps(resp.text)
    arr = jstring.split('length=')
    lstring = arr[1]
    arrl = lstring.split(';')
    print(lgene[i], arrl[0])






'''

fastadir = "D:/project/python/bioinfo/data/flybase/other"
with open(os.path.join(fastadir, "Ae2.json") , "w") as f:
    json.dump(data, f)
    '''
#save

'''
128up	FBgn0010339
Abp1	FBgn0036372
Acn	FBgn0263198
Acp98AB	FBgn0263597
Acsl	FBgn0263120
Act5C	FBgn0000042
Ahcy	FBgn0014455
alc	FBgn0260972
Alg1	FBgn0038552
Alh	FBgn0261238

fastadir = "D:/project/python/bioinfo/data/flybase/other"
with open(os.path.join(fastadir, "128up.json"), "r") as f:
    data = json.load(f)
#load
print("=========jsonObject.get(sequence)=============")
print(data['resultset']['result'][0]['sequence'])

fastadir = "D:/project/python/bioinfo/data/flybase/other"
bodyfiles = os.listdir(fastadir)
bodyseq={}
bodygene={}
idx = 0
'''
'''
for file in bodyfiles:
    filepath = os.path.join(fastadir, file)
    with open(filepath, "r") as f:
        data = json.load(f)
    bodygene[idx]= file.replace(".json","")
    bodyseq[idx] = data['resultset']['result'][0]['sequence']
    idx=idx+1
    print(data)

#print(data["contact"]["phone"]) # "123-456-7890"

print(bodyseq)
print(bodygene)

fastadir = "D:/project/python/bioinfo/data/flybase/npr"
nprfiles = os.listdir(fastadir)
nprseq={}
nprgene={}

idx = 0
for file in nprfiles:
    filepath = os.path.join(fastadir, file)
    seq_record = SeqIO.read(filepath, "fasta")
    nprgene[idx]= file.replace(".fasta","")
    nprseq[idx] = seq_record.seq
    idx=idx+1

print(nprseq)
print(nprgene)


aligner = Align.PairwiseAligner()

map ={}
score_size=25
s_size = 15
comp_size=30
score=0


for bi in bodyseq:
    print(bi)
    #print(bodyseq[bi])
    body_seq = bodyseq[bi]
    #print(body_seq)
    for ni in nprseq:
        print("=====npr=======",ni)
        #print(nprseq[ni])
        npr_seq =nprseq[ni]
        #print(npr_seq)
        #print(len(body_seq))
        #print(len(npr_seq))
        try:
            k = 0
            for i in range(int(len(body_seq) / s_size) + 1):
                for j in range(int(len(npr_seq) / s_size) + 1):
                    #print(body_seq[i * s_size:((i * s_size) + comp_size)])
                    #print(npr_seq[j * s_size:((j * s_size) + comp_size)])
                    alignments = aligner.align(body_seq[i * s_size:((i * s_size) + comp_size)],npr_seq[j * s_size:((j * s_size) + comp_size)])
                    #print("alignnem==",alignments)
                    score = aligner.score(body_seq[i * s_size:((i * s_size) + comp_size)],npr_seq[j * s_size:((j * s_size) + comp_size)])
                    #print(score)
                    if score >= score_size:
                        print(bodygene[bi], i * s_size, nprgene[ni], j * s_size, " score:", score, " count:", k)
                        map[bodygene[bi]+' '+nprgene[ni]] = (bodygene[bi]+":", i * s_size, nprgene[ni]+":", j * s_size, " score:", score, " count:", k)
                        print(alignments[0])
                        k = k + 1
                j = 0
        except:
            print("except:AstC_record,kr", score, k)


print("==========map===========")
print(map)
'''


