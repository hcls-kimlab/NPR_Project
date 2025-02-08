--create database
CREATE DATABASE drosphila DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci;

-- create user myuser/userpassord
create user 'myuser'@'localhost' identified by 'userpassord';

--data import to mysql

mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < gene.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < cell.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < motifregulon.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < motifregulonauc.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < tissue.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < tissue_cell.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < tissue_cell_gene.sql
mysql -umyuser -puserpassord --default-character-set=utf8 drosphila < tissue_gene.sql
