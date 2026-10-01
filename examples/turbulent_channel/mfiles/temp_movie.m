ANIMATE = ~false;

addpath '/Users/fringer/data/linux/suntans-trunk/mfiles';

datadir='../data';
suntansfile=[datadir,'/suntans.dat'];

nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
dt = getvalue(suntansfile,'dt');
nout = 1+nsteps/ntout;

figure(1);
clf;
for n=1:nout
    [x,z,d,T]=plotslice('T',datadir,n);
    [x,z,d,u]=plotslice('u',datadir,n);
    [x,z,d,w]=plotslice('w',datadir,n);        

    pcolor(x,z,T);
    colorbar;
    colormap jet;
    shading flat;
    axis image;
        
    drawnow;
end
