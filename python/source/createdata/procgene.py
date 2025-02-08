import numpy as np
import logging
import os
import sys
import loompy
import scanpy as sc
import numpy as np

print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))  # D:\data\slb\project\python
basedir = os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))  # D:\data\slb\project\python
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
print(os.path.dirname(os.path.abspath(os.path.dirname(os.path.dirname(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(__file__)))  # load_info

import sysinfo
import platform
import pymysql

datapath = sysinfo.datapath
loompath = sysinfo.loompath
print(loompath)

filelists = os.scandir(loompath)

conn = pymysql.connect(host=sysinfo.dburl, user=sysinfo.dbuser, password=sysinfo.dbpw,
                       db=sysinfo.database, charset='utf8')
curs = conn.cursor()


def getOid(objectname):
    print("getoid==========>")
    sql = " select objectvalue from oid where objectname='" + objectname + "'"
    curs.execute(sql)
    rows = curs.fetchall()
    row = rows[0]

    print("row==========>", row[0])

    sql = " update oid set objectvalue=" + str(int(row[0]) + 1) + " where objectname='" + objectname + "'"
    print(sql)
    curs.execute(sql)
    conn.commit()

    return row[0]


def isExistsGene(genename):
    sql = " select count(*) from gene where genename='" + genename + "'"
    curs.execute(sql)
    rows = curs.fetchall()
    row = rows[0]

    if (row[0] > 0):
        sql = " select geneid from gene where genename='" + genename + "'"
        curs.execute(sql)
        rows = curs.fetchall()
        row = rows[0]

        return row[0]
    else:
        return getOid("geneid")


def isExistsPart(partname):
    sql = " select count(*) from tissue where tissuename='" + partname + "'"
    curs.execute(sql)
    rows = curs.fetchall()
    print(rows)
    row = rows[0]

    print("row==>", row[0])

    if (int(row[0]) > 0):
        sql = " select tissueid from tissue where tissuename='" + partname + "'"
        curs.execute(sql)
        rows = curs.fetchall()
        row = rows[0]

        return row[0]
    else:
        return getOid("tissueid")


for file in filelists:
    filename = file.name
    print(filename)
    partname = os.path.splitext(filename)[0]
    print(partname)  # filename no extension
    fullpath = os.path.join(loompath, filename)
    pname = partname.replace("r_fca_biohub_", "")
    pname = pname.replace("s_fca_biohub_", "")
    pname = pname.replace("_10x", "")
    partid = isExistsPart(pname)
    print("partid[0]",partid)

    filedir = os.path.join(loompath, filename)
    print(filedir)
    ds = loompy.connect(filedir)
    genes = np.array(ds.row_attrs['Gene'])
    print(len(genes))
    genecount = len(genes)

    sql="update tissue set genecount="+str(genecount)+" where tissueid="+str(partid)
    curs.execute(sql)
    conn.commit()
    #part에 존재하는 gene의 갯수

    imgene = 0
    for val in genes:
        genename = val.decode('utf8')
        genename = genename.replace("'","´") #´'->디비에 저장이 안됨으로 변경후 저장
        print("genename",genename)
        geneid = isExistsGene(genename)

        print("geneid===>",geneid)

        sql = " insert into gene (geneid,genename) values ("+str(geneid)+",'"+genename+"')"
        sql = sql+" on duplicate key update tissuecount= tissuecount+1"

        print(sql)
        curs.execute(sql)
        conn.commit()

        #geneexpress



        sql = " insert into tissue_gene (tissueid, geneid) values ("+str(partid)+","+str(geneid)+")"
        sql = sql+" on duplicate key update icount= icount+1"


        print(sql)
        curs.execute(sql)
        conn.commit()
        imgene = imgene+1


ds.close()
conn.commit()
curs.close()
conn.close()

