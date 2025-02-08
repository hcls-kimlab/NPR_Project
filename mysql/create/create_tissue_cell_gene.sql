CREATE TABLE `tissue_cell_gene` (
  `tissueid` int(12) NOT NULL,
  `cellid` int(12) NOT NULL,
  `geneid` int(12) NOT NULL,
  `temp` varchar(200) DEFAULT NULL,
  `geneexpress` float(10,1) NOT NULL DEFAULT '1.0',
  `icount` int(12) NOT NULL DEFAULT '1',
  PRIMARY KEY (`tissueid`,`cellid`,`geneid`),
  KEY `index_cellid` (`cellid`),
  KEY `index_geneid` (`geneid`) USING BTREE,
  KEY `index_tissueid` (`tissueid`) USING BTREE
) ENGINE=InnoDB DEFAULT CHARSET=utf8