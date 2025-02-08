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
loompath="D:/project/python/bioinfo/data/loom1"
print(loompath)

filelists = os.scandir(loompath)

files=[]
for file in filelists:
    files.append(file.name)

print("sys.argv[1]==>",sys.argv[1])
ifile = int(sys.argv[1])

print("files===",files)
print("files[ifile]===",files[ifile])

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



#filename = file.name
filename = files[ifile]
print(filename)
partname = os.path.splitext(filename)[0]
print(partname)  # filename no extension
fullpath = os.path.join(loompath, filename)
pname = partname.replace("r_fca_biohub_", "")
pname = pname.replace("s_fca_biohub_", "")
pname = pname.replace("_10x", "")
partid = isExistsPart(pname)
print("partid[0]", partid)

filedir = os.path.join(loompath, filename)
print(filedir)
ds = loompy.connect(filedir)


#igenespcell = np.array(ds.col_attrs['n_genes'])
#print(len(igenespcell))
#genecount = len(genes)   #cell당 발견된 진의 수

cellidname = np.array(ds.col_attrs['CellID'])
print(len(cellidname))
cellcount = len(cellidname)  #cell 갯수


annotation = np.array(ds.col_attrs['annotation'])
print(len(annotation))
cellcount = len(annotation)  #cell 갯수

sql="update tissue set cellcount="+str(cellcount)+" where tissueid="+str(partid)
print(sql)
curs.execute(sql)
conn.commit()



#part에 존재하는 gene의 갯수

irow, icol = ds.shape
print("irow=>", irow)   #gene
print("icol=>", icol)   #cell
nparr = np.array(ds[:, :],dtype='f')

genes = np.array(ds.row_attrs['Gene'])
# print(len(genes))
genecount = len(genes)


icount=0
for cellid in cellidname :
    cellname = cellid.decode('utf8')
    #cellname = cellname.replace("'", "´")  # ´'->디비에 저장이 안됨으로 변경후 저장

#for icount in range(13423, len(cellidname)):

    #cellname = cellidname[icount].decode('utf8')
    cellname = cellname.replace("'","´") #´'->디비에 저장이 안됨으로 변경후 저장
    print("genename",cellname)

    cellid = isExistsCell(cellname)

    ann = annotation[icount].decode('utf8')
    ann = ann.replace("'", "´")  # ´'->디비에 저장이 안됨으로 변경후 저장

    print("cellid===>",cellid)

    sql = " insert into cell (cellid,cellname,annotation) values ("+str(cellid)+",'"+cellname+"','"+ann+"')"
    sql = sql+" on duplicate key update cellcount= cellcount+1"

    #print(sql)
    curs.execute(sql)
    conn.commit()

    sql = " insert into tissue_cell (tissueid, cellid) values ("+str(partid)+","+str(cellid)+")"
    sql = sql+" on duplicate key update icount= icount+1"

    #print(sql)
    curs.execute(sql)
    conn.commit()

    sql = " update cell set genecount= "+ str(ds.col_attrs['n_genes'][icount])
    sql = sql+" where cellid="+str(cellid)

    print(sql)
    curs.execute(sql)
    conn.commit()

    igcount = 0   #cell count

    for igcount in range(irow):
        #    #print(ds[j][2])

        if nparr[igcount, icount] >= 1.0:
            #print("nparr[igcount, icount]==>", nparr[igcount, icount])
            genename = genes[igcount].decode('utf8')
            genename = genename.replace("'", "´")  # ´'->디비에 저장이 안됨으로 변경후 저장
            geneid = isExistsGene(genename)
            geneexpress = nparr[igcount, icount]

            #print("geneexpress==>",geneexpress)

            sql = " insert into tissue_cell_gene (tissueid, cellid, geneid,geneexpress) values (" + str(
                partid) + "," + str(cellid) + "," + str(geneid) + "," + str(geneexpress) + ")"
            sql = sql + " on duplicate key update icount= icount+1"

            #print(sql)
            curs.execute(sql)
            conn.commit()

            igcount = igcount + 1



    icount=icount+1




ds.close()
conn.commit()
curs.close()
conn.close()

