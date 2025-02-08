CREATE TABLE `cell` (
  `cellid` int(12) NOT NULL COMMENT 'gene id',
  `cellname` varchar(300) DEFAULT NULL COMMENT 'gene name',
  `annotation` varchar(300) NOT NULL,
  `genecount` int(12) NOT NULL DEFAULT '1',
  `cellcount` int(12) NOT NULL DEFAULT '0',
  `sex` varchar(20) DEFAULT NULL,
  `age` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`cellid`) USING BTREE,
  KEY `index_name` (`cellname`),
  KEY `index_annotation` (`annotation`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8