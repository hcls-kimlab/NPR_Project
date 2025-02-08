---np

insert into prt_10_age_exp (tname, motifname, tfgenename, tfcount, npgenename, tfauc, 
npmotifregoccur, npmotifregweight, npgexpress, nprtype, age) 
select e.tname , np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  npgenename  , sum(np.motifauc) , sum(np.motifregoccur) , sum(np.motifregweight), sum(np.geneexpress),
 'np', np.age
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight , a.age
		from 
		(
			select cellid , age from np_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=17
		and a.cellid = i.cellid
		and b.tissueid = c.tissueid
		and c.reggeneid = d.geneid
		and d.geneid = i.geneid
		and c.tfgeneid = b.tfgeneid
		and d.used=1
		and d.genetype='p'
		and b.tissueid= i.tissueid
		and d.geneid= e.geneid
) np ,
cell a ,
gene b , 
gene c ,
tissue e 
where np.cellid = a.cellid
and np.tfgeneid = b.geneid
and np.npgeneid = c.geneid
and np.tissueid = e.tissueid
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid , np.age
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid , np.age

--npr

insert into prt_10_age_exp (tname, motifname, tfgenename, tfcount, nprgenename, tfauc, 
nprmotifregoccur, nprmotifregweight, nprgexpress, nprtype, age) 
select e.tname , np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  nprgenename  ,  sum(np.motifauc) ,  sum(np.motifregoccur) , sum(np.motifregweight) 
,sum(np.geneexpress), 'npr', np.age
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight , a.age
		from 
		(
			select cellid , age from npr_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=17
		and a.cellid = i.cellid
		and b.tissueid = c.tissueid
		and c.reggeneid = d.geneid
		and d.geneid = i.geneid
		and c.tfgeneid = b.tfgeneid
		and d.used=1
		and d.genetype='r'
		and b.tissueid= i.tissueid
		and d.geneid= e.geneid
) np ,
cell a ,
gene b , 
gene c ,
tissue e 
where np.cellid = a.cellid
and np.tfgeneid = b.geneid
and np.npgeneid = c.geneid
and np.tissueid = e.tissueid
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid , np.age
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid , np.age


--np 존재

select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount ,  nprtype, age
from
(
select a.* from  prt_10_age_exp a
where a.nprtype='np' and a.age = '3'
) aa
group by aa.motifname , aa.npgenename , aa.nprtype , aa.age

--npr 존재
select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype, age
from
(
select a.* from  prt_10_age_exp a
where a.nprtype='npr' and a.age = '3'
     
) aa
group by aa.motifname , aa.nprgenename , aa.nprtype , aa.age



select aa.tname , replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount ,  nprtype, age
from
(
select a.* from  prt_10_age_exp a
where a.nprtype='np' and a.age = '3'
) aa
group by aa.tname, aa.motifname , aa.npgenename , aa.nprtype , aa.age

--npr 존재
select aa.tname, replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype, age
from
(
select a.* from  prt_10_age_exp a
where a.nprtype='npr' and a.age = '3'
     
) aa
group by aa.tname, aa.motifname , aa.nprgenename , aa.nprtype , aa.age



select npp.tname, npp.motifname motifname, npp.tfcount , nprr.tfcount, npp.genename, npp.npcount, nprr.genename , nprr.nprcount from
(

	select tname , motifname, sum(tfcount) tfcount , npgenename genename, sum(npmotifregoccur) npcount , pid
	from 
	(
		select a.tname , a.motifname , a.tfcount, a.npgenename , a.npmotifregoccur , c.pid from prt_10_exp a ,
		gene b , (
		select a.geneid , b.pid , a.genetype , a.used from nprtgene a,pairprgene b
		where a.used = 1
		and a.used = b.used
		and a.geneid = b.geneid
		and a.genetype='p'
		) c
		where a.nprtype='np'
		and a.npgenename = b.genename
		and b.geneid = c.geneid
		and c.genetype='p'
		and c.used=1
		and a.tname in ('body','head')

	) np
	group by tname , motifname , npgenename , npmotifregoccur , pid
	order by tname , pid
) npp left join
(
	select tname , motifname, sum(tfcount) tfcount, nprgenename genename, sum(nprmotifregoccur) nprcount , pid
	from 
	(
	select a.tname , a.motifname , a.tfcount, a.nprgenename , a.nprmotifregoccur , c.pid from prt_10_exp a ,
	gene b , (
	select a.geneid , b.pid , a.genetype , a.used from nprtgene a,pairprgene b
	where a.used = 1
	and a.used = b.used
	and a.geneid = b.geneid
	and a.genetype='r'
	) c
	where a.nprtype='npr'
	and a.nprgenename = b.genename
	and b.geneid = c.geneid
	and c.genetype='r'
	and c.used=1
	and a.tname in ('body','head')
	) npr
	group by tname , motifname , pid
    order by tname , motifname, pid
)
nprr on (nprr.tname = npp.tname
and nprr.motifname = npp.motifname
and npp.pid = nprr.pid
)


select npp.tname, npp.motifname motifname, npp.genename, npp.npcount, nprr.genename , nprr.nprcount from
(
	select tname , motifname, sum(tfcount) tfcount, nprgenename genename, sum(nprmotifregoccur) nprcount , pid
	from 
	(
	select a.tname , a.motifname , a.tfcount, a.nprgenename , a.nprmotifregoccur , c.pid from prt_10_exp a ,
	gene b , (
	select a.geneid , b.pid , a.genetype , a.used from nprtgene a,pairprgene b
	where a.used = 1
	and a.used = b.used
	and a.geneid = b.geneid
	and a.genetype='r'
	) c
	where a.nprtype='npr'
	and a.nprgenename = b.genename
	and b.geneid = c.geneid
	and c.genetype='r'
	and c.used=1
	and a.tname in ('body','head')
	) npr
	group by tname , motifname , pid
    order by tname , pid , motifname
    

	
) nprr right join
(
  
    select tname , motifname, sum(tfcount), npgenename genename, sum(npmotifregoccur) npcount , pid
	from 
	(
		select a.tname , a.motifname , a.tfcount, a.npgenename , a.npmotifregoccur , c.pid from prt_10_exp a ,
		gene b , (
		select a.geneid , b.pid , a.genetype , a.used from nprtgene a,pairprgene b
		where a.used = 1
		and a.used = b.used
		and a.geneid = b.geneid
		and a.genetype='p'
		) c
		where a.nprtype='np'
		and a.npgenename = b.genename
		and b.geneid = c.geneid
		and c.genetype='p'
		and c.used=1
		and a.tname in ('body','head')

	) np
	group by tname , motifname , npgenename , npmotifregoccur , pid
	order by tname , pid ,motifname
)
npp on 
(nprr.tname = npp.tname
and nprr.motifname = npp.motifname
and npp.pid = nprr.pid
)

