%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: profplot.m
% Description: Plot the data produced by the profile output in suntans.
%
% Oliver Fringer
% Stanford University
% 07/15/06
%
% This mfile will plot the data output along the transect specified 
% by the parameters in suntans.dat.  An example of the parameters
% in suntans.dat is given below:
%
% ProfileVariables	default          # Output u, s, s0, and h
% DataLocations	        dataxy.dat       # dataxy.dat contains column x-y data
% ProfileDataFile	profdata.dat     # Information about profiles is in profdata.dat
% ntoutProfs		10               # Output profile data every 10 time steps
% NkmaxProfs		100              # Only output the top 10 z-levels
% numInterpPoints	1                # Output data at the three nearest neighbors to each input point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataX, dataY, u, v, rtime] = profvels(datadir, n)

dirname = datadir;
EMPTY = 999999;

fname = [dirname,'/profdata.dat'];

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

dataX = dataXY(1:2:end);
dataY = dataXY(2:2:end);

dzh = 0.5*(dz(1:end-1)+dz(2:end));
z = zeros(Nkmax,1);
z(1) = -dz(1)/2;
z(2:Nkmax) = z(1) - cumsum(dzh);
dZ = dz*ones(1,length(dataX));

x = sqrt((dataX-dataX(1)).^2+(dataY-dataY(1)).^2)';
X = ones(Nkmax,1)*x;
Z = z*ones(1,numTotalDataPoints);

nout=nsteps/ntoutProfs;

iplot=fix(numTotalDataPoints/2);
kplot=fix(Nkmax/2);
% 
% %%% added by Brooke
% %%% also %'d out every s0
% 
% iplot = find(dataX>90 & dataX<90.5 & dataY>-30 & dataY<-29); %point in embayment
% %iplot = find(dataX>(90+150) & dataX<(91+150) & dataY>-40 & dataY<-39.5); %point in embayment
% if kplot == 0
%     kplot=1;
% end

fname = [dirname,'/u.dat.prof'];
ufid = fopen(fname,'rb');

rtime = n*ntoutProfs*dt;
for n=1:n
    udata = reshape(fread(ufid,3*Nkmax*numInterpPoints*numTotalDataPoints,'float64'),...
              Nkmax,numInterpPoints,numTotalDataPoints,3);
    udata(find(udata==EMPTY))=nan;


  
    u = squeeze(udata(:,1,:,1));
    v = squeeze(udata(:,1,:,2));

    uplot(n) = u(iplot, kplot);
    vplot(n) = v(iplot, kplot);
end

fclose(ufid);

% dT = ntoutProfs*dt;
% t = dT*[1:nout];

% subplot(4,1,1)
% plot(t,hplot);
% set(gca,'xticklabel','');
% ylabel('h(t)');

% figure(2)
% 
subplot(4,1,2)
plot(t/3600,uplot);
set(gca,'xticklabel','');
ylabel('u(t)');
% 
% subplot(4,1,3)
% plot(t/3600,sqrt(vplot.^2+uplot.^2));
% set(gca,'xticklabel','');
% ylabel('vel mag(t)');
% 
% subplot(4,1,4)
% plot(t/3600,vplot);
% xlabel('t');
% ylabel('v(t)');

end



  



