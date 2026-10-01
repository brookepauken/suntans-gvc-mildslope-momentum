LOAD = ~false;
%datadir = '../data_highres';
%datadir = '../data';
EMPTY = 999999;

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';
datadir = '../data_turbtest';

fname = [datadir,'/profdata.dat'];

suntansfile=[datadir,'/suntans.dat'];

profvels('data_turbtest', 1)