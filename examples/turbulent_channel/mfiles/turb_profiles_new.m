% addpath '/Users/fringer/data/linux/suntans-trunk/mfiles';
% 
% datadir = '../data_highres';

%Findr_Nx_minx(10, 32, .1)

clear all; close all;

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';
datadir = '../data_turbtest';

cd '/Users/brookepauken/suntans-gvc-mild-slope-2025/examples/turbulent-channel/mfiles';

datadir = '~/new_adv_tests/copy_newadvect_3dcoarse_p0ab3'


EMPTY = 999999;

fname = [datadir,'/profdata.dat'];

suntansfile=[datadir,'/suntans.dat'];

fid = fopen(fname,'rb');
numTotalDataPoints = fread(fid,1,'int32');
numInterpPoints = fread(fid,1,'int32');
Nkmax = fread(fid,1,'int32');
nsteps = fread(fid,1,'int32');
ntoutProfs = fread(fid,1,'int32');  
dt = fread(fid,1,'float64');
dz = fread(fid,Nkmax,'float64');
dataIndices = fread(fid,numTotalDataPoints,'int32');
dataXY = fread(fid,2*numTotalDataPoints,'float64');
xv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
yv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
fclose(fid);

H = sum(dz);
kappa = 0.4;
ustar = 0.02;
nu = getvalue(suntansfile,'nu');

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
num_avg = 3;

nout=nsteps/ntoutProfs;

fname = [datadir,'/u.dat.prof'];
ufid = fopen(fname,'rb');
fname = [datadir,'/s.dat.prof'];

uread = fread(ufid,'float64');
fclose(ufid);
Nread = length(uread);
Nperstep = 3*Nkmax*numInterpPoints*numTotalDataPoints;
nout = floor(length(uread)/Nperstep);
uread = uread(1:nout*Nperstep);

udata = reshape(uread,Nkmax,numInterpPoints,numTotalDataPoints,3,nout);

uall = squeeze(udata(:,1,:,1,:));
vall = squeeze(udata(:,1,:,2,:));
wall = squeeze(udata(:,1,:,3,:));

ubar = zeros(Nkmax,1);
vbar = zeros(Nkmax,1);
wbar = zeros(Nkmax,1);
nstart = max(1,nout - floor(num_avg*n_tau));
nfinish = nout;
 
num_avg_shift = 1;
nstart = max(1,nout - floor(num_avg_shift*n_tau) - floor(num_avg*n_tau));
nfinish = min(nout,nout - floor(num_avg_shift*n_tau));

for n=nstart:nfinish
    for p=1:numTotalDataPoints
        ubar = ubar + uall(:,p,n);
        vbar = vbar + vall(:,p,n);
        wbar = wbar + wall(:,p,n);
    end
end
num_total = (1+nout-nstart)*numTotalDataPoints;
ubar = ubar/num_total;
vbar = vbar/num_total;
wbar = wbar/num_total;

uprime = zeros(Nkmax,1);
vprime = zeros(Nkmax,1);
wprime = zeros(Nkmax,1);
uw = zeros(Nkmax,1);
for n=nstart:nfinish
    for p=1:numTotalDataPoints
        uprime = uprime + (uall(:,p,n)-ubar).^2;
        vprime = vprime + (vall(:,p,n)-vbar).^2;
        wprime = wprime + (wall(:,p,n)-wbar).^2;
        uw = uw + (uall(:,p,n)-ubar).*(wall(:,p,n)-wbar);        
    end
end
uprime = sqrt(uprime/num_total);
vprime = sqrt(vprime/num_total);
wprime = sqrt(wprime/num_total);
uw = uw/num_total;

fontsize=14;

figure(1)
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
       'location','northeast');
xlabel('u_{rms}/u_*');
ylabel('z^+');
set(gca,'fontsize',14);
sgtitle('Averaged over 3 periods')

U_th = 1/H*sum(dz.*u_th);
kplot = 10;
prof_index = 2;
u_series = squeeze(uall(kplot,prof_index,:))/U_th;
v_series = squeeze(vall(kplot,prof_index,:))/U_th;
w_series = squeeze(wall(kplot,prof_index,:))/U_th;
ubar_series = sum(squeeze(uall(:,prof_index,:)).*(dz*ones(1,nout)))/H/U_th;
t = [1:nout]*dt*ntoutProfs;
figure(2)
plot(t/tau,u_series,'k-',t/tau,v_series,'r-',t/tau,w_series,'b-',t/tau,ubar_series,'m-');
xlabel('$u_* t/H$');
ylabel('$u$');
set(gca,'fontsize',14);
grid on 
hold on 
plot(t(nstart)/tau*ones(1, 100), linspace(-0.4, 1.6, 100), 'k--', 'LineWidth', 1.5)
plot(t(nfinish)/tau*ones(1, 100), linspace(-0.4, 1.6, 100), 'k--', 'LineWidth', 1.5)
xlim([0, t(nout)/tau]+1)
legend('$u_{inst}$/mean($u_{log}$)','$v_{inst}$/mean($u_{log}$)','$v_{inst}$/mean($u_{log}$)','ubar/mean($u_{log}$)',...
       '','','location','best');
ylim([0.6, 1.6])