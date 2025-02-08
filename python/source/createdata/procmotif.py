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
def isExistsCell(cellname):
    sql = " select count(*) from cell where cellname='" + cellname + "'"
    curs.execute(sql)
    rows = curs.fetchall()
    row = rows[0]

    if (row[0] > 0):
        sql = " select cellid from cell where cellname='" + cellname + "'"
        curs.execute(sql)
        rows = curs.fetchall()
        row = rows[0]

        return row[0]
    else:
        return getOid("cellid")


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


sql = " select geneid , genename from gene order by geneid "
curs.execute(sql)
rows = curs.fetchall()
dictgene = {row[1]: row[0] for row in rows}

ifile = 0
for file in filelists:
    filename = file.name
    print(filename)
    partname = os.path.splitext(filename)[0]
    print(partname)  # filename no extension
    pname = partname.replace("r_fca_biohub_", "")
    pname = pname.replace("s_fca_biohub_", "")
    pname = pname.replace("_10x", "")

    tissueid = isExistsPart(pname)
    print("tissueid===>",tissueid)
    fullpath = os.path.join(loompath, filename)


    filedir = os.path.join(loompath, filename)
    print(filedir)
    ds = loompy.connect(filedir)

    tf = [x for x, y in ds.row_attrs['MotifRegulons'].dtype.fields.items()]
    print(len(tf))
    print(tf[0])  # motif full name
    print(tf[0].index('_('))
    print(tf[0][:tf[0].index('_(')])  # motif gene name
    # 'MotifRegulonGeneWeights
    print(ds.row_attrs['MotifRegulons'][tf[0]])  # matrix를 가져온다..

    genes = np.array(ds.row_attrs['Gene'])

    for i in range(len(tf)):
        print(tf[i])
        motifname = tf[i]
        motifgenename = tf[i][:tf[i].index('_(')]  # motif gene name
        motifgenename = motifgenename.replace("'", "´")  # ´'->디비에 저장이 안됨으로 변경후 저장
        try:
            motifgeneid = dictgene[motifgenename]
        except:
            continue

        motifregeulons = np.transpose(np.nonzero(ds.row_attrs['MotifRegulons'][motifname]))

        for regrow in motifregeulons:
            print(regrow[0])  #gene matrix
            reggenename = genes[regrow[0]].decode('utf8')
            reggenename = reggenename.replace("'", "´")  # ´'->디비에 저장이 안됨으로 변경후 저장

            print(reggenename)
            try:
                reggeneid = dictgene[reggenename]
            except:
                reggeneid=0
            print(reggenename,reggeneid)

            motifcount=0
            motifweight=0.0
            motifoccur = 0

            try:
                motifcount= ds.row_attrs['MotifRegulons'][motifname][regrow[0]]
            except:
                motifcount = 0

            try:
                motifweight = ds.row_attrs['MotifRegulonGeneWeights'][motifname][regrow[0]]
            except:
                motifweight = 0.0

            try:
                motifoccur = ds.row_attrs['MotifRegulonGeneOccurrences'][motifname][regrow[0]]
            except:
                motifoccur = motifcount


            sql = " insert into motifregulon (tissueid,motifname,tfgeneid,reggeneid,motifreg,motifregweight,motifregoccur) "
            sql = sql+" values ("+str(tissueid)+",'"+str(motifname)+"',"+str(motifgeneid)+","
            sql = sql+str(reggeneid)+","+str(motifcount)+","+str(motifweight)+","+str(motifoccur)+")"
            sql = sql+" on duplicate key update nprttype='d'"

            print(motifcount)
            print(motifweight)  # matrix를 가져온다..
            print(motifoccur)  # matrix를 가져온다..

            print(sql)
            curs.execute(sql)
            conn.commit()


ds.close()
conn.commit()
curs.close()
conn.close()

