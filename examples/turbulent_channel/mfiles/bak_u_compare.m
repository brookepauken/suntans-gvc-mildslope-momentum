clear all; close all;

MOVIE=~false;

% addpath '/Users/fringer/data/linux/suntans-trunk/mfiles';
% datadir='../data';

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';

datadir = '../data_turbtest';


suntansfile=[datadir,'/suntans.dat'];

nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
Nkmax = getvalue(suntansfile,'Nkmax');
nout = nsteps/ntout;

H = 10;
kappa = 0.4;
ustar = 0.02;
nu = getvalue(suntansfile,'nu');

fname = [datadir,'/profdata.dat'];
fid = fopen(fname,'rb');
numTotalDataPoints = fread(fid,1,'int32');
numInterpPoints = fread(fid,1,'int32');
nsteps = fread(fid,1,'int32');
ntoutProfs = fread(fid,1,'int32');  
dt = fread(fid,1,'float64');
dz = fread(fid,Nkmax,'float64');

dzh = 0.5*(dz(1:end-1)+dz(2:end));
z = zeros(Nkmax,1);
z(1) = -dz(1)/2;
z(2:Nkmax) = z(1) - cumsum(dzh); % u is locate half-way between w faces
zw = dz(1)-cumsum(dz(1:Nkmax)); % w is located at faces

z0 = nu/9/ustar;
u_th = ustar/kappa*log((z+H)/z0);
z_plus = ustar*(H+z)/nu;
u_th(z_plus<11) = ustar^2*(z(z_plus<11)+H)/nu;
tau = H/ustar;
n_tau = floor(tau/dt/ntoutProfs);
num_avg = 2;

u0 = 0.2831;

%nout=1;


ubar = zeros(Nkmax,1);
vbar = zeros(Nkmax,1);
wbar = zeros(Nkmax,1);


for n=1:nout
    fprintf('On %d of %d\n',n,nout);
    
    [x,z,dv,u]=plotslice('u',datadir,n);
    [x,z,dv,v]=plotslice('v',datadir,n);    
    [x,z,dv,w]=plotslice('w',datadir,n);
    [~,~, ~, h]=plotslice('h',datadir,n); 
    
    ubar = ubar + mean(u,2);
    vbar = vbar + mean(v,2);
    wbar = wbar + mean(w,2);
 
    if(MOVIE)
        figure(1)
        totvel = sqrt(u.^2 + v.^2 + w.^2)./u0;
        pcolor(x, z, totvel)
        colorbar;
        a = colorbar;
        clim([prctile(totvel, 1, "all"), prctile(totvel, 95, "all")]);
        a.Label.String = '$u_{mag}/u0$';
        a.Label.Interpreter = 'latex';
        a.Label.FontSize=16;
        xlabel('x (m)')
        ylabel('z (m)')
        shading flat
        daspect([1 1 1])

        % figure(1)
        % plot(x, h);
        % xlabel('x (m)')
        % ylabel('h (m)')


        % figure(2)
        % pcolor(x, z, w)
        % colorbar;
        % xlabel('x (m)')
        % ylabel('z (m)')
        % shading interp
        % %pause;
    end

end

ubar = ubar/nout;
vbar = vbar/nout;
wbar = wbar/nout;

uprime = zeros(Nkmax,1);
vprime = zeros(Nkmax,1);
wprime = zeros(Nkmax,1);

%%
for n=1:nout
    fprintf('On %d of %d\n',n,nout);
    
    [x,z,dv,u]=plotslice('u',datadir,n);
    [x,z,dv,v]=plotslice('v',datadir,n);    
    [x,z,dv,w]=plotslice('w',datadir,n);
    
    Nc = length(ubar(:))/Nkmax;
    up = u - ubar*ones(1,Nc);
    vp = v - vbar*ones(1,Nc);
    wp = w - wbar*ones(1,Nc);        

    uprime = uprime + up.^2;
    vprime = vprime + vp.^2;
    wprime = wprime + wp.^2;


end
uprime = sqrt(uprime/nout);
vprime = sqrt(vprime/nout);
wprime = sqrt(wprime/nout);

figure(2)
subplot(1,2,1)
semilogx(z_plus,ubar/ustar,'ko',...
         z_plus,u_th/ustar,'k-');
xlabel('z^+');
ylabel('u^+');
set(gca,'fontsize',14);

subplot(1,2,2)
plot(uprime/ustar,z_plus,'k-',vprime/ustar,z_plus,'r-',wprime/ustar,z_plus,'b-',...
     uw/ustar^2,z_plus,'m-',[-1 0],[0 max(z_plus)],'k--');
legend('u_{rms}','v_{rms}','w_{rms}','Reynolds stress',...
       'location','northwest');
xlabel('u_{rms}/u_*');
ylabel('z^+');
set(gca,'fontsize',14);
