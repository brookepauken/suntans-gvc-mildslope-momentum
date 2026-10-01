LOAD = ~false;
%datadir = '../data_highres';
datadir = '../data';
EMPTY = 999999;

fname = [datadir,'/profdata.dat'];

suntansfile=[datadir,'/suntans.dat'];

H = 10;
kappa = 0.4;
ustar = 0.02;
nu = getvalue(suntansfile,'nu');
Nkmax = getvalue(suntansfile,'Nkmax');
dt = getvalue(suntansfile,'dt');
nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');

dz = load([datadir,'/vertspace.dat']);
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
n_tau = floor(tau/dt/ntout);
num_avg = 4;

c=load([datadir,'/cells.dat']);
Nc0=length(c(:,1));
xv=c(:,2);
yv=c(:,3);

file_info = dir([datadir,'/u.dat']);
filesize = file_info.bytes;
nout=filesize/8/Nkmax/Nc0/3;

W = 25;
Ny = 50;
dy = W/Ny;
j = Ny/2;
yc = (j-1/2)*dy;
yc = 0.5*W;
Nx=200;
j=30;
is = [1:Nc0];%Nx*(j-1)+[1:Nx];find(abs(yv-yc)<dy/2);
Nc = length(is)

if(LOAD)
    ubar = zeros(Nkmax,1);
    vbar = zeros(Nkmax,1);
    wbar = zeros(Nkmax,1);
    nstart = max(1,nout - floor(num_avg*n_tau));
    fid = fopen([datadir,'/u.dat'],'rb');
    for n=nstart:nout
        fprintf('On %d of %d\n',n,nout);

        %        fseek(fid,Nkmax*3*Nc0*(n-1)*8,'bof');
        %        uvw = reshape(fread(fid,Nkmax*3*Nc0,'float64'),Nc0,Nkmax,3);
        uvw = reshape(getcdata(fid,Nkmax*3*Nc0,n,'float64'),Nc0,Nkmax,3);
        %        [x,z,dv,uvw]=plotslice('U',datadir,n);
        %        [x,z,dv,v]=plotslice('v',datadir,n);    
        %        [x,z,dv,w]=plotslice('w',datadir,n);
        u = squeeze(uvw(is,:,1))';
        v = squeeze(uvw(is,:,2))';
        w = squeeze(uvw(is,:,3))';
        
        ubar = ubar + u;
        vbar = vbar + v;
        wbar = wbar + w;
    end
    num_total = (1+nout-nstart);
    ubar = mean(ubar,2)/num_total;
    vbar = mean(vbar,2)/num_total;
    wbar = mean(wbar,2)/num_total;
    
    uprime = zeros(Nkmax,1);
    vprime = zeros(Nkmax,1);
    wprime = zeros(Nkmax,1);
    uw = zeros(Nkmax,1);
    for n=nstart:nout
        fprintf('On %d of %d\n',n,nout);

        fseek(fid,Nkmax*3*Nc0*(n-1)*8,'bof');
        uvw = reshape(fread(fid,Nkmax*3*Nc0,'float64'),Nc0,Nkmax,3);
        %        [x,z,dv,uvw]=plotslice('U',datadir,n);
        %        [x,z,dv,v]=plotslice('v',datadir,n);    
        %        [x,z,dv,w]=plotslice('w',datadir,n);
        u = squeeze(uvw(is,:,1))';
        v = squeeze(uvw(is,:,2))';
        w = squeeze(uvw(is,:,3))';
        
        %        Nc = length(ubar(:))/Nkmax;
        up = u - ubar*ones(1,Nc);
        vp = v - vbar*ones(1,Nc);
        wp = w - wbar*ones(1,Nc);        
        
        uprime = uprime + up.^2;
        vprime = vprime + vp.^2;
        wprime = wprime + wp.^2;
        uw = uw + up.*wp;
    end
    fclose(fid);
    
    uprime = sqrt(mean(uprime,2)/num_total);
    vprime = sqrt(mean(vprime,2)/num_total);
    wprime = sqrt(mean(wprime,2)/num_total);
    uw = mean(uw,2)/num_total;
end

fontsize=14;

figure(1)
subplot(1,2,1)
semilogx(z_plus,ubar/ustar,'r-',...
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


