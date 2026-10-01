addpath '/Users/fringer/data/linux/suntans-trunk/mfiles';

datadir = '../data3';
EMPTY = 999999;

suntansfile=[datadir,'/suntans.dat'];

dz = load([datadir,'/vertspace.dat']);

H = sum(dz);
kappa = 0.4;
ustar = 0.02;
nu = getvalue(suntansfile,'nu');
Nkmax = getvalue(suntansfile,'Nkmax');

dzh = 0.5*(dz(1:end-1)+dz(2:end));
z = zeros(Nkmax,1);
z(1) = -dz(1)/2;
z(2:Nkmax) = z(1) - cumsum(dzh); % u is locate half-way between w faces

z0 = nu/9/ustar;
u = ustar/kappa*log((z+H)/z0);
ubar = 1/H*sum(dz.*u);

fprintf('Ubar = %.4f\n',ubar);
