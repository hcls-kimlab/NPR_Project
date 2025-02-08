--np,npr,np+npr +tf each tissue expression
--

--np_npr annotation 10

select e.tname , a.annotation,  np_npr.npmotifname,np_npr.nprmotifname, b.genename  tfgenename , count(np_npr.tfgeneid) tfcount , 
c.genename  npgenename  , count(np_npr.npgeneid) npcount , d.genename  nprgenename , count(np_npr.nprgeneid) nprcount , np_npr.npauc , np_npr.nprauc ,
avg(f.motifregoccur) , avg(f.motifregweight), avg(g.motifregoccur) , avg(g.motifregweight)  from 
(
select np.tissueid, np.cellid , np.tfgeneid , np.motifname npmotifname , np.geneid npgeneid , npr.motifname nprmotifname , npr.geneid nprgeneid , np.motifauc npauc , npr.motifauc nprauc from
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,  d.geneid , e.pid ,  b.motifauc 
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
	) a,
	motifregulonauc b ,
	motifregulon c ,
	nprtgene d ,
    pairprgene e,
	tissue_cell_gene i
	where a.cellid = b.cellid
    and b.tissueid=8
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
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,   d.geneid ,e.pid , b.motifauc
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
	) a,
	motifregulonauc b ,
	motifregulon c ,
	nprtgene d ,
    pairprgene e,
	tissue_cell_gene i
	where a.cellid = b.cellid
	and a.cellid = i.cellid
    and b.tissueid=8
	and b.tissueid = c.tissueid
	and c.reggeneid = d.geneid
	and d.geneid = i.geneid
	and c.tfgeneid = b.tfgeneid
	and d.used=1
	and d.genetype='r'
	and b.tissueid= i.tissueid
    and d.geneid = e.geneid
    ) npr    
where np.tissueid = npr.tissueid
and np.cellid = npr.cellid
and np.tfgeneid = npr.tfgeneid
and np.pid = npr.pid
) np_npr ,
cell a ,
gene b , 
gene c ,
gene d ,
tissue e ,
motifregulon f ,
motifregulon g
where np_npr.cellid = a.cellid
and f.tissueid=8
and g.tissueid=8
and np_npr.tfgeneid = b.geneid
and np_npr.tfgeneid = f.tfgeneid
and np_npr.npmotifname = f.motifname
and np_npr.tfgeneid = g.tfgeneid
and np_npr.nprmotifname = g.motifname
and np_npr.npgeneid = c.geneid
and np_npr.npgeneid = f.reggeneid
and np_npr.nprgeneid = d.geneid
and np_npr.nprgeneid = g.reggeneid
and np_npr.tissueid = e.tissueid
group by e.tname , a.annotation, np_npr.tfgeneid , np_npr.npgeneid
order by e.dispord , a.annotation , np_npr.tfgeneid , np_npr.npgeneid


--np+npr annotation 10 
select e.tname , a.annotation,  np_npr.npmotifname, np_npr.nprmotifname, b.genename  tfgenename , count(np_npr.tfgeneid) tfcount , 
c.genename  npgenename  , count(np_npr.npgeneid) npcount , d.genename  nprgenename , count(np_npr.nprgeneid) nprcount , np_npr.npauc , np_npr.nprauc ,
avg(f.motifregoccur) , avg(f.motifregweight), avg(npgexpress) ,  avg(g.motifregoccur) , avg(g.motifregweight) , avg(nprgexpress) from 
(
select np.tissueid, np.cellid , np.tfgeneid , np.motifname npmotifname , np.geneid npgeneid , np.geneexpress npgexpress, npr.geneexpress nprgexpress,  npr.motifname nprmotifname , npr.geneid nprgeneid , np.motifauc npauc , npr.motifauc nprauc from
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,  d.geneid , e.pid ,  b.motifauc , i.geneexpress
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
	) a,
	motifregulonauc b ,
	motifregulon c ,
	nprtgene d ,
    pairprgene e,
	tissue_cell_gene i
	where a.cellid = b.cellid
    and b.tissueid=3
	and a.cellid = i.cellid
	and b.tissueid = c.tissueid
	and c.reggeneid = d.geneid
	and d.geneid = i.geneid
	and c.tfgeneid = b.tfgeneid
	and c.motifname = b.motifname
	and d.used=1
	and d.genetype='p'
	and b.tissueid= i.tissueid
    and d.geneid= e.geneid
) np ,
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,   d.geneid ,e.pid , b.motifauc , i.geneexpress
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
	) a,
	motifregulonauc b ,
	motifregulon c ,
	nprtgene d ,
    pairprgene e,
	tissue_cell_gene i
	where a.cellid = b.cellid
	and a.cellid = i.cellid
    and b.tissueid=3
	and b.tissueid = c.tissueid
	and c.reggeneid = d.geneid
	and d.geneid = i.geneid
	and c.tfgeneid = b.tfgeneid
    and c.motifname = b.motifname
	and d.used=1
	and d.genetype='r'
	and b.tissueid= i.tissueid
    and d.geneid = e.geneid
    ) npr    
where np.tissueid = npr.tissueid
and np.cellid = npr.cellid
and np.tfgeneid = npr.tfgeneid
and np.pid = npr.pid
) np_npr ,
cell a ,
gene b , 
gene c ,
gene d ,
tissue e ,
motifregulon f ,
motifregulon g
where np_npr.cellid = a.cellid
and f.tissueid=3
and g.tissueid=3
and np_npr.tfgeneid = b.geneid
and np_npr.tfgeneid = f.tfgeneid
and np_npr.npmotifname = f.motifname
and np_npr.nprmotifname = g.motifname
and np_npr.tfgeneid = g.tfgeneid
and np_npr.npgeneid = c.geneid
and np_npr.npgeneid = f.reggeneid
and np_npr.nprgeneid = d.geneid
and np_npr.nprgeneid = g.reggeneid
and np_npr.tissueid = e.tissueid
group by e.tname , a.annotation, np_npr.npmotifname, np_npr.nprmotifname , np_npr.tfgeneid , np_npr.npgeneid
order by e.dispord , a.annotation , np_npr.npmotifname, np_npr.nprmotifname, np_npr.tfgeneid ,  np_npr.tfgeneid , np_npr.npgeneid
 

select * from motifregulon

--np 10 
select e.tname , a.annotation,  np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  npgenename  , count(np.npgeneid) npcount , avg(np.motifauc) , avg(np.geneexpress),
avg(np.motifregoccur) , avg(np.motifregweight) 
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight
		from 
		(
			select cellid from np_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=3
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
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid

1,2,4,5,6,7,9,10,11,12,13,14,15,16,17
insert into 
--npr 10 
select e.tname , a.annotation,  np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  npgenename  , count(np.npgeneid) npcount , avg(np.motifauc) , avg(np.geneexpress),
avg(np.motifregoccur) , avg(np.motifregweight) 
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight
		from 
		(
			select cellid from npr_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=8
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
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid





--np+npr tissue 10 

select a.* from prt_10_exp a , tissue b
where a.tname = b.tname
order by b.dispord , a.nprtype, a.motifname , a.npgenename , a.nprgenename

select * from tissue

--np+npr+tf 데이타...

insert into prt_10_exp (tname, motifname, tfgenename, tfcount, npgenename, nprgenename, tfauc, 
npmotifregoccur, npmotifregweight, npgexpress, nprmotifregoccur, nprmotifregweight, nprgexpress, nprtype) 
select e.tname ,  np_npr.nprmotifname, b.genename  tfgenename , count(np_npr.tfgeneid) tfcount , 
c.genename  npgenename  , d.genename  nprgenename , np_npr.npauc ,
CAST(avg(f.motifregoccur)AS signed integer)  , avg(f.motifregweight),avg(npgexpress) ,  
CAST(avg(g.motifregoccur)AS signed integer)  , avg(g.motifregweight) , avg(nprgexpress) ,'np_npr' from 
(
select np.tissueid, np.cellid , np.tfgeneid , np.motifname npmotifname , np.geneid npgeneid , np.geneexpress npgexpress, npr.geneexpress nprgexpress,  npr.motifname nprmotifname , npr.geneid nprgeneid , np.motifauc npauc , npr.motifauc nprauc from
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,  d.geneid , e.pid ,  b.motifauc , i.geneexpress
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
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
(
	select b.tissueid , a.cellid , b.tfgeneid, b.motifname,   d.geneid ,e.pid , b.motifauc , i.geneexpress
	from 
	(
		select j.cellid from npr_10_cell  j,
		np_10_cell k
		where j.cellid = k.cellid
	) a,
	motifregulonauc b ,
	motifregulon c ,
	nprtgene d ,
    pairprgene e,
	tissue_cell_gene i
	where a.cellid = b.cellid
	and a.cellid = i.cellid
    and b.tissueid=17
	and b.tissueid = c.tissueid
	and c.reggeneid = d.geneid
	and d.geneid = i.geneid
	and c.tfgeneid = b.tfgeneid
	and d.used=1
	and d.genetype='r'
	and b.tissueid= i.tissueid
    and d.geneid = e.geneid
    ) npr    
where np.tissueid = npr.tissueid
and np.cellid = npr.cellid
and np.tfgeneid = npr.tfgeneid
and np.pid = npr.pid
) np_npr ,
cell a ,
gene b , 
gene c ,
gene d ,
tissue e ,
motifregulon f ,
motifregulon g
where np_npr.cellid = a.cellid
and f.tissueid=17
and g.tissueid=17
and np_npr.tfgeneid = b.geneid
and np_npr.tfgeneid = f.tfgeneid
and np_npr.npmotifname = f.motifname
and np_npr.nprmotifname = g.motifname
and np_npr.tfgeneid = g.tfgeneid
and np_npr.npgeneid = c.geneid
and np_npr.npgeneid = f.reggeneid
and np_npr.nprgeneid = d.geneid
and np_npr.nprgeneid = g.reggeneid
and np_npr.tissueid = e.tissueid
group by e.tname , e.tname , np_npr.npmotifname, np_npr.nprmotifname , np_npr.tfgeneid , np_npr.npgeneid
order by e.dispord , e.tname , np_npr.npmotifname, np_npr.nprmotifname, np_npr.tfgeneid ,  np_npr.tfgeneid , np_npr.npgeneid
 










select * from motifregulon

select * from 
--np tissue 10 
select * from prt_10_exp where nprtype='np'

insert into prt_10_exp (tname, motifname, tfgenename, tfcount, 
npgenename, tfauc, npgexpress, npmotifregoccur, npmotifregweight, nprtype) 
select e.tname ,  np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  npgenename  , avg(np.motifauc) , avg(np.geneexpress), CAST(avg(np.motifregoccur) AS signed integer) , avg(np.motifregweight) , 'np'
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight
		from 
		(
			select cellid from np_10_cell
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
group by e.tname , e.tname , np.motifname, np.tfgeneid , np.npgeneid
order by e.dispord , e.tname , np.motifname, np.tfgeneid , np.npgeneid



--npr tissue 10 
select * from prt_10_exp
select * from tissue
alter table prt_10_exp_back modify tname varchar(100);

insert into prt_10_exp (tname, motifname, tfgenename, tfcount, npgenename, nprgenename, tfauc, npmotifregoccur, npmotifregweight, npgexpress, nprmotifregoccur, 
nprmotifregweight, nprgexpress, nprtype) 
select e.tname ,  np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , '' ,
c.genename  nprgenename  , avg(np.motifauc) tfauc ,0,0,0, CAST(avg(np.motifregoccur) AS signed integer) nprmotifregoccur ,  
avg(np.motifregweight) nprmotifregweight , avg(np.geneexpress) nprgeneexpress , 'npr'
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight
		from 
		(
			select cellid from npr_10_cell
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
group by e.tname , e.tname, np.motifname, np.tfgeneid , np.npgeneid
order by e.dispord , e.tname , np.motifname, np.tfgeneid , np.npgeneid

select * from npr

select a.tname, motifname , sum(tfcount) from prt_10_exp  a 
, tissue b
where nprtype='np' and
a.tname = b.tname
group by b.dispord , motifname 


select a.tname, motifname , sum(tfcount) from prt_10_exp  a 
, tissue b
where nprtype='npr' and
a.tname = b.tname
group by b.dispord , motifname 


select a.tname , npgenename , sum(npmotifregoccur) from prt_10_exp  a 
, tissue b
where nprtype='np' and 
a.tname = b.tname
group by b.dispord , npgenename 

select a.tname , nprgenename , sum(nprmotifregoccur) from prt_10_exp  a 
, tissue b
where nprtype='npr' and 
a.tname = b.tname
group by  b.dispord ,nprgenename 

--whole tissue
select a.* from prt_10_exp a, tissue c
where  a.nprtype='np'
and a.tname = c.tname
and not exists
( 
   select 1 from prt_10_exp b
    where nprtype='np_npr'
    and b.motifname=a.motifname
)
order by c.dispord , a.motifname, npgenename

select a.* from prt_10_exp a , tissue c
where  a.nprtype='npr'
and a.tname = c.tname
and not exists
( 
	select 1 from prt_10_exp b
    where nprtype='np_npr'
    and b.motifname=a.motifname
)
order by c.dispord , a.motifname, nprgenename



--all

    select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount , nprtype
    from 
    (
		select a.* from prt_10_exp a, tissue c
		where  a.nprtype='np'
		and a.tname = c.tname
		and not exists
		( 
		   select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		)
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.npgenename , aa.nprtype
    
    
	select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype
    from 
    (
		select a.* from prt_10_exp a, tissue c
		where  a.nprtype='npr'
		and a.tname = c.tname
		and not exists
		( 
		   select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		)
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.nprgenename , aa.nprtype