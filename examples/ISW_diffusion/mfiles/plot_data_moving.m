%

% Brooke Pauken
%

close all; clear all;

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';

ANIMATE=true;
SAVEMOVIE=true;

cwave=1.7; %or whatever it is


%set datadir name and string to append to start of movie files
datadir = '../data_diffusion_test';
modifier = 'example-movie';


shortrun= 1; %1 for quick runs where we want timestep as title, 0 otherwise
plotvar = 1; %1 for salinity, 2 for temp, 3 for velocity, 4 for free surface, 5 for pressure, 7 for contours of s,
                %8 for voriticty, 9 for Ri, 11 for nut
 

suntansfile=[datadir,'/suntans.dat'];

nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
z0 = getvalue(suntansfile,'z0B');
CdBprop = getvalue(suntansfile,'CdB');
nu = getvalue(suntansfile,'nu'); %viscosity
dt = getvalue(suntansfile,'dt');
beta = getvalue(suntansfile,'beta');
nout = 1+floor(nsteps/ntout);

coordsys = getvalue(suntansfile,'vertcoord');
if(coordsys==2 || coordsys==4)
    iso =1;
else
    iso = 0;
end

%cwave = getvalue(suntansfile,'cwave');

time_vec = 1:ntout:nsteps+1; 
time_vec = time_vec*dt;

figure(1)
clf;



% Load the grid to determine the length of the domain and the cell centers
p = load([datadir,'/points.dat']);
xp = p(:,1);
L = max(xp);
cells = load([datadir,'/cells.dat']);
%edges = load([datadir,'/edges.dat']);
xv = cells(:,2);
dx = xv(2) - xv(1); %find dx by finding distance between two xs
yv = cells(:,3);
dy = yv(1)*2; %dy just a constant, since y is at center of cell we need this x 2
dz = load([datadir,'/vertspace.dat']);
dzb = dz(end);
nk = length(dz);
nx = length(xv);

W = dy; %width is constant 
D = sum(dz);


if(SAVEMOVIE && ANIMATE)
    if plotvar==1
        vidtit = ['./results/splot_', modifier];
    elseif plotvar==2
        vidtit = ['./results/Tplot_', modifier];
    elseif plotvar==3
        vidtit = ['./results/uplot_', modifier];
    elseif plotvar==4
        vidtit = ['./results/fs_', modifier];
    elseif plotvar==7
        vidtit = ['./results/splotconts_', modifier];
    elseif plotvar==8
        vidtit = ['./results/vorplot_', modifier];
    elseif plotvar==9
        vidtit = ['./results/richardsonplot_', modifier];
    elseif plotvar==11
        vidtit = ['./results/nutplot_', modifier];       
    else
        vidtit = ['/results/otherplot_', modifier];
    end
    vid=VideoWriter(vidtit, 'MPEG-4');
    vid.FrameRate = 30;
    vid.Quality = 90;
    open(vid)

end



% Stratification is given by a tanh distribution
D = round(sum(dz));         % Depth
rho0 = 1000;         % Reference density
g = 9.81;            % Gravity
delta = 0.15*D; %0.1*D;       % Pycnocline thickness


ubot_max = zeros(1, nout);
u_time_av = zeros(nk, nx);
w_time_av = zeros(nk, nx);

cnt = 0;

ssum = zeros(1, ntout);
eta_mid=zeros(1, nout);
max_ubot=zeros(1, nout);
L_sep = zeros(1, nout);
L_main = zeros(1, nout);
L_front = zeros(1, nout);


u_nearbot = zeros(1, nout);
w_nearbot = zeros(1, nout);
z_B_offline = zeros(1, length(xv));


umid_nearbot=zeros(1, nout);
umid_neartop=zeros(1, nout);

w_nearbed=zeros(1, nout);
u_nearbed=zeros(1, nout);
x_spot=zeros(1, nout);


for n=1:1:nout
    cnt = cnt+1;
    %fprintf('On %d of %d\n',n,nout);
    [x,z, dv, s]=plotslice_moving('s',datadir,n);   
    [~,~, ~, T]=plotslice_moving('T',datadir,n);   
    [~,~, ~, u]=plotslice_moving('u',datadir,n); 
    [~,~, ~, w]=plotslice_moving('w',datadir,n); 
    [~,~, ~, h]=plotslice_moving('h',datadir,n); 
    [~, ~, ~, q]=plotslice_moving('q', datadir,n);
    if(plotvar==11)
        [~, ~, ~, nut]=plotslice_moving('nut', datadir,n);
    end

  
    u_time_av = u_time_av+u;
    w_time_av = w_time_av+w;


    %get dzzs
    fid = fopen([datadir,'/vert_dzz.dat'],'rb');
    fseek(fid,8*nk*nx*(n-1),'bof');
    dzz = reshape(fread(fid,nk*nx,'float64'),nx,nk);
    fclose(fid);

    depthdata = load([datadir,'/depth.dat-voro'],'-ascii');
    dv = depthdata(:,3);

    zB = ((dv+h) - sum(dzz, 2));


    %if you want the pcolor plots to span the whole depth
    if(iso)
        z(1, :) = 0;
        z(end, :) = -dv;
    end

    ubot_max(n) = max(abs(u(end, :)));
    l2s = sqrt(sum(u.^2, 1));
    

    if n == 1

        conts = contourc(x(1, :), z(:, 1), s, 1);
        mean_z = mean(conts(2, 2:end));
        %salin_thresh = conts(1, 1);
        salin_thresh = 0;
        u_base = u;
        u_old = u;
        w_base = w;
        s_base = s;
        s_old = s;

        z_min = min(min(z));
        z_max = max(max(z(end, :)));
    

        if(ANIMATE)
            figure(1)
            %fontsize=18;
            clf;
            [curlz, cav]=curl(x, z, u-mean(mean(u)), w-mean(mean(w)));
            zplot = 0 - D; 
           
            if(plotvar==1) %salinty
               stmp = (s)*beta*rho0+rho0; %convert s to density
             
               pcolor(x/1000,z,stmp);

               colorbar;
               a=colorbar;
               a.Label.String = '$\rho$ ($kg/m^3$)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
               pbaspect([6 1 1])

               set(gcf, 'Position', [140   521   921   259]);

               ylabel('z (m)')
               xlabel('x (km)')
               shading flat;
            elseif(plotvar==11) %nu_T
               pcolor(x/1000,z,nut);
               colorbar;
               a=colorbar;
               a.Label.String = '$\nu_T$ ($kg/m^3$)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
            elseif(plotvar==2) %temperature
               pcolor(x/1000,z, T);
               colorbar;
               a=colorbar;
               a.Label.String = 'T';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               shading flat;
           elseif(plotvar==3) %velocity
               pcolor(x/1000, z, u);
               shading flat
               a=colorbar;
               a.Label.String = '$u$ (m/s)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
               shading flat
               cmocean('balance', 'pivot',0)
               fontsize(16, 'points');
           elseif(plotvar==5) %pressure things
               ptot = rho0*g*(h'+z);

               p0 = rho0*g*(h'+z); %base 
               pb = g*cumtrapz(rho(z).*dzz');
                
               rhoprime = s*rho0+rho0 - rho(z) - rho0;

               pprime = g*cumtrapz(rhoprime.*dzz');

               ptot = p0+pb+pprime + q;
               pplot = pprime + q;
               %pplot = q;

               pcolor(x/1000,z,pplot);

               colorbar;
               a=colorbar;
               a.Label.String = 'P';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               shading interp;

    
               [dpdx, dphdx] = pgradx(x, z, s, h, dzz, g, rho0, q);
    
               subplot(2, 1, 1)
               pcolor(x/1000,z/D,dpdx);
               colorbar;
               a=colorbar;
               a.Label.String = 'dP/dx';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylim([-1, -.9])
               shading interp
               colormap jet
               pbaspect([8 1 1])
               ylabel('z/H')
               %clim([-15, 15])
               ax = gca;
               ax.FontSize = 16; 
    
               hold on
               subplot(2, 1, 2)
               pcolor(x/1000,z/D,curlz);
               colorbar;
               a=colorbar;
               a.Label.String = 'vorticity';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylim([-1, -.9])   
               clim([-.2, .2])
               shading interp
               colormap jet
               pbaspect([8 1 1])
               ylabel('z/H')
               xlabel('x (km)')
               ax = gca;
               ax.FontSize = 16; 
            elseif(plotvar==7) %contours of s
                contour(x/1000,z, s, 30, 'k');
                ylabel('z (m)', 'FontSize', 16);
                xlabel('x (km)', 'FontSize', 16);
            elseif(plotvar==8) %near-bed vorticity
                pcolor((x-5000)/D,z/D,curlz);
                colorbar;
                ylim([-1, -.9]);         
                clim([-.2, .2]);
                a=colorbar;
                a.Label.String = 'vorticity ($1/s$)';
                colormap turbo;
                shading interp;
                xlabel('x/H');
                ylabel('z/H');
                pbaspect([10 1 1])
                hold on 
                cmocean('balance', 'pivot',0)
            elseif(plotvar==9) %richardson number
                % Density gradient and buoyancy frequency N
                drho_dz = zeros(nk ,nx);
                dudz = zeros(nk, nx);
                
                z_ri = zeros(nk, nx);
                dzz=dzz';
                z_ri(1, :) = -dzz(1, :)/2;
                for k=2:nk
                    z_ri(k, :) = z_ri(k-1, :) - 0.5*(dzz(k-1, :)+dzz(k, :));
                end
                
                % higher_modes.m needs the grid to increase with increasing index, unlike the
                % suntans grid, which is the opposite
                %z_ri = z_ri(end:-1:1, :);

                for k=2:nk-1
                    drho_dz(k, :) = ((s(k+1, :)*1000+1000)-(s(k-1, :)*1000+1000))./(z_ri(k+1, :)-z_ri(k-1, :));
                    %drho_dz(k, :) = ((s(k+1, :))-(s(k-1, :)))/(z(k+1, :)-z(k-1, :));
                    %drho_dz(k, :) = (rho_bar(k+1)-rho_bar(k-1))/(z(k+1)-z(k-1));

                    dudz(k, :) = ((u(k+1, :))-(u(k-1, :)))./(z_ri(k+1, :)-z_ri(k-1, :));

                    drho_dz(k, :) = ((s(k+1, :)*1000+1000)-(s(k, :)*1000+1000))./(z_ri(k+1, :)-z_ri(k, :));
                    dudz(k, :) = ((u(k+1, :))-(u(k, :)))./(z_ri(k+1, :)-z_ri(k, :));

        
                end
                drho_dz(1, :) = 1.5*drho_dz(2, :)-0.5*drho_dz(3, :);
                drho_dz(nk, :) = 1.5*drho_dz(nk-1, :)-0.5*drho_dz(nk-2, :);

                dudz(1, :) = 1.5*dudz(2, :)-0.5*dudz(3, :);
                dudz(nk, :) = 1.5*dudz(nk-1, :)-0.5*dudz(nk-2, :);

                N2 = (-g/rho0).*drho_dz;
                Ri = N2./(dudz.^2);

                pcolor(x/1000,z, Ri)
                shading flat
                colorbar
                a=colorbar;
                a.Label.String = 'Richardson Number';
                colormap jet;
                shading interp;
                xlabel('x (km)');
                ylabel('z');
                clim([0, 1]);
            elseif(plotvar==4) %free surface
               plot(x(1, :)/1000,h*1000);
               ylabel('Free surface height (mm)', 'FontSize', 16);
               xlabel('x (km)', 'FontSize', 16);
            end
    
            if(plotvar~=4 && plotvar~=10 && plotvar~=12)
                % inds_above = find(zplot>(-D+dz(end)/2));
                % zplot_above = zplot;
                % zplot_above(~inds_above)=NaN;
                % zplot_below=zplot;
                % zplot_below(inds_above)=NaN;
                % 
                % x_above=x(1, :);
                % x_above(~inds_above)=NaN;
                % x_below=x(1, :);
                % x_below(inds_above)=NaN;
                % 
                % plot(x_above/1000, zplot_above, 'w')
                % plot(x_below/1000, zplot_below, 'k')
                if(shortrun==0)
                    titlestr = strcat('$tc/L$ =  ', num2str(round((time_vec(n))/(L/cwave), 2)));
                elseif(shortrun==1)
                    titlestr = strcat('$t$ =  ', num2str(round((time_vec(n)/3600), 2)));  
                    titlestr = strcat(titlestr, ' hours'); 
                end

                if(plotvar~=5 && plotvar~=12)
                    title(titlestr, 'FontSize', 16);
                else
                    subplot(2, 1, 1)
                    title(titlestr, 'FontSize', 16);
                end
                
              
            end

            set(gcf, 'Color', 'white');
            if(plotvar==3 || plotvar==8)
                figitem = figure(1);
                figitem.Position = [300 300 1000 300];
            elseif(plotvar==5)
                figitem = figure(1);
                figitem.Position =[1 275 1512 587];
            end
            clear fontsize;
            fontsize(16, 'points');   
            drawnow;
            pause;
            if(SAVEMOVIE && plotvar~=6)
                F = getframe(figure(1)); 
                writeVideo(vid, F);

            end
        end
        
    end


    [rows, cols] = size(s);
    zvec = z(:, 1)';
    
    if n==1
        running_sum = u;
    else
        running_sum = running_sum + u;
    end
    
    av_u = running_sum/n;
    u_old = u;
    ssum(n) = sum(sum(s));
    eta_mid(n)=D+z(end, length(xv)/2);

    if(ANIMATE)
        figure(1)
        clf;
        [curlz, cav]=curl(x, z, u, w);

        if(plotvar==1) %salinty
               stmp = (s)*beta*rho0+rho0; %convert s to density
             
               pcolor(x/1000,z,stmp);

               colorbar;
               a=colorbar;
               a.Label.String = '$\rho$ ($kg/m^3$)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
               pbaspect([6 1 1])

               set(gcf, 'Position', [140   521   921   259]);

               ylabel('z (m)')
               xlabel('x (km)')
               shading flat;
            elseif(plotvar==11) %nu_T
               pcolor(x/1000,z,nut);
               colorbar;
               a=colorbar;
               a.Label.String = '$\nu_T$ ($kg/m^3$)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
            elseif(plotvar==2) %temperature
               pcolor(x/1000,z, T);
               colorbar;
               a=colorbar;
               a.Label.String = 'T';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               shading flat;
           elseif(plotvar==3) %velocity
               pcolor(x/1000, z, u);
               shading flat
               a=colorbar;
               a.Label.String = '$u$ (m/s)';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylabel('z (m)')
               xlabel('x (km)')
               shading flat
               cmocean('balance', 'pivot',0)
               fontsize(16, 'points');
           elseif(plotvar==5) %pressure things
               ptot = rho0*g*(h'+z);

               p0 = rho0*g*(h'+z); %base 
               pb = g*cumtrapz(rho(z).*dzz');
                
               rhoprime = s*rho0+rho0 - rho(z) - rho0;

               pprime = g*cumtrapz(rhoprime.*dzz');

               ptot = p0+pb+pprime + q;
               pplot = pprime + q;
               %pplot = q;

               pcolor(x/1000,z,pplot);

               colorbar;
               a=colorbar;
               a.Label.String = 'P';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               shading interp;

    
               [dpdx, dphdx] = pgradx(x, z, s, h, dzz, g, rho0, q);
    
               subplot(2, 1, 1)
               pcolor(x/1000,z/D,dpdx);
               colorbar;
               a=colorbar;
               a.Label.String = 'dP/dx';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylim([-1, -.9])
               shading interp
               colormap jet
               pbaspect([8 1 1])
               ylabel('z/H')
               %clim([-15, 15])
               ax = gca;
               ax.FontSize = 16; 
    
               hold on
               subplot(2, 1, 2)
               pcolor(x/1000,z/D,curlz);
               colorbar;
               a=colorbar;
               a.Label.String = 'vorticity';
               a.Label.Interpreter = 'latex';
               a.Label.FontSize=16;
               ylim([-1, -.9])   
               clim([-.2, .2])
               shading interp
               colormap jet
               pbaspect([8 1 1])
               ylabel('z/H')
               xlabel('x (km)')
               ax = gca;
               ax.FontSize = 16; 
            elseif(plotvar==7) %contours of s
                contour(x/1000,z, s, 30, 'k');
                ylabel('z (m)', 'FontSize', 16);
                xlabel('x (km)', 'FontSize', 16);
            elseif(plotvar==8) %near-bed vorticity
                pcolor((x-5000)/D,z/D,curlz);
                colorbar;
                ylim([-1, -.9]);         
                clim([-.2, .2]);
                a=colorbar;
                a.Label.String = 'vorticity ($1/s$)';
                colormap turbo;
                shading interp;
                xlabel('x/H');
                ylabel('z/H');
                pbaspect([10 1 1])
                hold on 
                cmocean('balance', 'pivot',0)
            elseif(plotvar==9) %richardson number
                % Density gradient and buoyancy frequency N
                drho_dz = zeros(nk ,nx);
                dudz = zeros(nk, nx);
                
                z_ri = zeros(nk, nx);
                dzz=dzz';
                z_ri(1, :) = -dzz(1, :)/2;
                for k=2:nk
                    z_ri(k, :) = z_ri(k-1, :) - 0.5*(dzz(k-1, :)+dzz(k, :));
                end
                
                % higher_modes.m needs the grid to increase with increasing index, unlike the
                % suntans grid, which is the opposite
                %z_ri = z_ri(end:-1:1, :);

                for k=2:nk-1
                    drho_dz(k, :) = ((s(k+1, :)*1000+1000)-(s(k-1, :)*1000+1000))./(z_ri(k+1, :)-z_ri(k-1, :));
                    %drho_dz(k, :) = ((s(k+1, :))-(s(k-1, :)))/(z(k+1, :)-z(k-1, :));
                    %drho_dz(k, :) = (rho_bar(k+1)-rho_bar(k-1))/(z(k+1)-z(k-1));

                    dudz(k, :) = ((u(k+1, :))-(u(k-1, :)))./(z_ri(k+1, :)-z_ri(k-1, :));

                    drho_dz(k, :) = ((s(k+1, :)*1000+1000)-(s(k, :)*1000+1000))./(z_ri(k+1, :)-z_ri(k, :));
                    dudz(k, :) = ((u(k+1, :))-(u(k, :)))./(z_ri(k+1, :)-z_ri(k, :));

        
                end
                drho_dz(1, :) = 1.5*drho_dz(2, :)-0.5*drho_dz(3, :);
                drho_dz(nk, :) = 1.5*drho_dz(nk-1, :)-0.5*drho_dz(nk-2, :);

                dudz(1, :) = 1.5*dudz(2, :)-0.5*dudz(3, :);
                dudz(nk, :) = 1.5*dudz(nk-1, :)-0.5*dudz(nk-2, :);

                N2 = (-g/rho0).*drho_dz;
                Ri = N2./(dudz.^2);

                pcolor(x/1000,z, Ri)
                shading flat
                colorbar
                a=colorbar;
                a.Label.String = 'Richardson Number';
                colormap jet;
                shading interp;
                xlabel('x (km)');
                ylabel('z');
                clim([0, 1]);
            elseif(plotvar==4) %free surface
               plot(x(1, :)/1000,h*1000);
               ylabel('Free surface height (mm)', 'FontSize', 16);
               xlabel('x (km)', 'FontSize', 16);
            end

        
        if(shortrun==0)
            titlestr = strcat('$tc/L$ =  ', num2str(round((time_vec(n))/(L/cwave), 2)));
        elseif(shortrun==1)
            titlestr = strcat('$t$ =  ', num2str(round((time_vec(n)/3600), 2)));  
            titlestr = strcat(titlestr, ' hours');
        end
        if(plotvar~=5 && plotvar~=12 && plotvar~=10)
            title(titlestr, 'FontSize', 16);
        end
     
        
        if(plotvar~=4 & plotvar~=12 & plotvar~=10)
            a.Label.Interpreter = 'latex';
            a.Label.FontSize=16;
            shading interp;
        end

        
        shading flat;
        drawnow;
        clear fontsize;
        fontsize(16, 'points');   
        if(SAVEMOVIE) 
            F = getframe(figure(1)); 
            writeVideo(vid, F);
        end

    end

    s_old = s;
    s_RMSE(n) = sqrt(sum(sum((s-s_base).^2)));

end


if(SAVEMOVIE)
    close(vid)
end

