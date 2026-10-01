MOVIE=~false;
addpath '/Users/fringer/data/linux/suntans-trunk/mfiles';

datadir='../data_y';
suntansfile=[datadir,'/suntans.dat'];

nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
nout = nsteps/ntout;

xmax = 100;

figure(1)
clf;
subplot(4,1,1);
hold on;
subplot(4,1,2);
hold on;
subplot(4,1,3);
hold on;
subplot(4,1,4);
hold on;

shade_type = 'flat';

if(MOVIE)
    for n=1:nout
        fprintf('On %d of %d\n',n,nout);

        [x,z,dv,u]=plotslice_y('u',datadir,n);
        [x,z,dv,w]=plotslice_y('w',datadir,n);        
        [x,z,dv,T]=plotslice_y('T',datadir,n);        

        subplot(4,1,1)
        cla;
        pcolor(x,z,u);
        colorbar;
        colormap jet;
        shading(shade_type);
        xlabel('x (m)');
        ylabel('z (m)');
        title('u velocity (m/s)');

        subplot(4,1,2)
        cla;
        pcolor(x,z,w);
        colorbar;
        colormap jet;
        shading(shade_type);
        xlabel('x (m)');
        ylabel('z (m)');
        title('w velocity (m/s)');        

        subplot(4,1,3)
        cla;
        quiver(x,z,u,w,'k-');
        plot(x(1,:),-dv,'k-');
        xlabel('x (m)');
        ylabel('z (m)');
        title('Quiver');
        
        subplot(4,1,4)
        cla;
        pcolor(x,z,T);        
        colorbar;
        colormap jet;
        shading(shade_type);
        %        plot([0 xmax],[-7.5 -7.5],'r--',...
        %             [0 xmax],[-12.5 -12.5],'r--');         
        xlabel('x (m)');
        ylabel('z (m)');
        title('Scalar');
        
        drawnow;
    end
else
    [x,z,d,u]=plotslice('u',datadir,nout);
end

xv = x(1,:);
xplot = [25,50,75];
iplot = zeros(size(xplot));
for m=1:length(xplot)
    iplot(m) = min(find(abs(xv-xplot(m))==min(abs(xv-xplot(m)))));
end

u0 = 0.1;
figure(2)
plot(u(:,iplot)/u0,z(:,iplot),'k-');
xlabel('u/u_0');
ylabel('z (m)');
