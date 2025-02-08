CREATE TABLE `motifregulonauc` (
  `tissueid` int NOT NULL DEFAULT '0',
  `cellid` int NOT NULL DEFAULT '0',
  `motifname` varchar(100) NOT NULL,
  `tfgeneid` int NOT NULL DEFAULT '0',
  `motifauc` float(10,9) NOT NULL DEFAULT '0.000000000',
  `nprttype` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`tissueid`,`cellid`,`tfgeneid`,`motifname`),
  KEY `index_motifname` (`motifname`),
  KEY `index_tfgeneid` (`tfgeneid`),
  KEY `index_cellid` (`cellid`),
  KEY `index_tissueid` (`tissueid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3