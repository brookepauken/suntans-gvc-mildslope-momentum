%

% Brooke Pauken
%


%close all; clear all;

%exampledir = '/Users/brookepauken/suntans-gvc-mild-slope/examples/nudge_DJL_wave/mfiles/';
%exampledir = '/Users/brookepauken/suntans-gvc-mild-slope/examples/periodic_DJL_wave/sherlock_out/';


%cd(exampledir);
ANIMATE=true;
ENERGYPLOT = ~true;
SAVEMOVIE=true;
SAVEUBOT=true; %make true to save figures and bottom velocities, will plot either way

addpath '/home/groups/fringer/bpauken/mfiles';

if ~exist('sh_val','var')
    sh_val = 2;
    nk = 200;
end
nx = 800;
dt = 0.9121;
sh_val = 2;

%set up datadir name and naming convention for figures
multiplier = sh_val/10;

%modifier = ['pproject_h1_', num2str(multiplier*10)];
%modifier = ['srun_Nx_', num2str(nx)];
%modifier = ['dt_', num2str(dt)];
%modifier = ['dt_moving_dt_', num2str(dt), '_copy'];
%week='092524/';
%modifier='copy_data_movingbed';
%modifier = ['copy_bl_test4_correct', num2str(nx)];
%modifier = 'copy_data_stretch_nude_20_r1.005';
%datadir=['../data_', modifier];
%datadir=['/Users/brookepauken/suntans-gvc-mild-slope/examples/periodic_DJL_wave/sherlock_out/data_', modifier];
%datadir=['/Users/brookepauken/suntans-gvc-mild-slope/examples/nudge_DJL_wave/', modifier];
%datadir=['/Volumes/ExHard/sherlock_out/', week, modifier];
%datadir=['/Volumes/ExHard/sherlock_out/', modifier];
%datadir=['/scratch/users/bpauken/suntans_output_temp/', modifier];

suntansfile=[datadir,'/suntans.dat'];


nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
dt = getvalue(suntansfile,'dt');
nout = 1+floor(nsteps/ntout);

cwave = getvalue(suntansfile,'cwave');
%cwave = 0; %comment out for wave frame

time_vec = 0:ntout:nsteps;
time_vec = time_vec*dt;

figure(1)
clf;

% Load the grid to determine the length of the domain and the cell centers
p = load([datadir,'/points.dat']);
xp = p(:,1);
L = max(xp);
cells = load([datadir,'/cells.dat']);
xv = cells(:,2);
dx = xv(2) - xv(1); %find dx by finding distance between two xs
yv = cells(:,3);
dy = yv(1)*2; %dy just a constant, since y is at center of cell we need this x 2
dz = load([datadir,'/vertspace.dat']);
dzb = dz(end);

W = dy; %width is constant 


if(SAVEMOVIE && ANIMATE)
    vidtit = [datadir, '/results/splot_', modifier];
    vid=VideoWriter(vidtit, 'MPEG-4');
    vid.FrameRate = 30;
    vid.Quality = 90;
    open(vid)

    vidtit = [datadir, '/results/uplot_', modifier];
    vid2=VideoWriter(vidtit, 'MPEG-4');
    vid2.FrameRate = 30;
    vid2.Quality = 90;
    open(vid2)

    vidtit = [datadir, '/results/vorplot_', modifier];
    vid3=VideoWriter(vidtit, 'MPEG-4');
    vid3.FrameRate = 30;
    vid3.Quality = 90;
    open(vid3)
end

% vid=VideoWriter(vidtit, 'MPEG-4');
% vid.FrameRate = 30;
% vid.Quality = 90;
% open(vid)


g = 9.8;
rho0 = 1000;

ubot_max = zeros(1, nout);
etta_max = zeros(1, nout);
amps_mid = zeros(1, nout);
amps_max = zeros(1, nout);
error = zeros(1, nout);
w_error = error;
max_error = zeros(1, nout);
av_u_plot = zeros(1, nout);
av_u_diff = zeros(1, nout);

width = zeros(1, nout);

Eb = zeros(1, ntout);
Ep = zeros(1, ntout);
Ke = zeros(1, ntout);

uflip = NaN(1, nout);

cnt = 0;

ssum = zeros(1, ntout);


for n=1:nout
    cnt = cnt+1;
    %fprintf('On %d of %d\n',n,nout);
    %[x,z, dv, s]=plotslice_moving('s',datadir,n);   
    %[~,~, ~, u]=plotslice_moving('u',datadir,n); 
    %[~,~, ~, w]=plotslice_moving('w',datadir,n); 

    [x,z, s]=plotslice('s',datadir,n);   
    [~,~, u]=plotslice('u',datadir,n); 
    [~,~, w]=plotslice('w',datadir,n);

    ubot_max(n) = max(abs(u(end, :)));

    if n == 1
        D = abs(z(end, 1));
        mean_z = multiplier*D; %height of pycnocline for calculating amplitudes
        conts = contourc(x(1, :), z(:, 1), s, 1);
        mean_z = mean(conts(2, 2:50));
        salin_thresh = conts(1, 1);
        u_base = u;
        u_old = u;
        w_base = w;
        s_base = s;
        
    end

    [rows, cols] = size(s);
    zvec = z(:, 1)';

    amps_max_xs = zeros(1, cols);

    for ival = 1:cols
        zval = interp1(s(:, ival), zvec, salin_thresh);
        amps_max_xs(ival) = zval;
    end

    zmean = max(amps_max_xs); %logically this should maybe be h1 instead...
    % figure;
    % plot(xv, amps_max_xs)

    % figure(1)
    % plot(x(1, :), amps_max_xs)
    % hold off

    if sh_val < 5
        zmean = max(amps_max_xs);
    else
        zmean = min(amps_max_xs);
    end

    amps_max_xs = amps_max_xs - zmean;
    % figure;
    % plot(xv, amps_max_xs)

    % 
    % amps = conts(2, 2:end) - mean_z;
    % amps_mid(n) = interp1(conts(1, 2:end), amps, 2500);
    % amps_max(n) = max(abs(amps)); 

    amps_max(n) = max(abs(amps_max_xs));

    width(n) = (1/amps_max(n))*trapz(amps_max_xs)*dx; %characteristic width

    xspot=241;
    %uxspot=u(:, xspot)+cwave;
    uxspot=u(:, xspot);
    
    indneg=find(uxspot<=0);
    uxspot=uxspot(indneg(end):end);
    indreverse = find(uxspot>0);
    if isempty(indreverse)
        uflip(n)=0;
    else
        uflip(n)=D+zvec(indreverse(1)+indneg(end)-1);
    end

    % if uflip(n)>(-D*.55)
    %     uflip(n)=0;
    % end


    if n==1
        running_sum = u;
    else
        running_sum = running_sum + u;
    end
    
    av_u = running_sum/n;
    av_u_plot(n) = u(20, 100) - av_u(20, 100);
    av_u_plot(n) = u(20, 100) - u_base(20, 100);
    u_diff = u - u_old;
    av_u_diff(n) = mean(mean(u_diff));
    u_old = u;

    
    error(n) = sqrt(sum(sum((u-u_base).^2))/numel(u)); %try without abs
    w_error(n) = sqrt(sum(sum((w-w_base).^2))/numel(w)); %try without abs
    max_error(n) = max(max(abs(u-u_base)));
    

    ssum(n) = sum(sum(s));


    

    if(ANIMATE)
        % if(n>1)
        %    delete(botplot);
        % end
        figure(1)
        %[curlz, cav]=curl(x, z, u+cwave, w);
        [curlz, cav]=curl(x, z, u, w);
        pcolor(x/D,z/D,s);
        %ylim([-300, -275]);
        hold on
        %botplot= plot(x(1,:)/1000,-dv,'k-');
        titlestr = strcat('time = ', num2str(time_vec(n)*cwave/D), ' t*c/H');
        title(titlestr);
        colorbar;
        a=colorbar;
        a.Label.String = 'salinity';
        colormap jet;
        shading flat;
        xlabel('x/H');
        ylabel('z/H');
        
        drawnow;
        if(SAVEMOVIE)
            F = getframe(figure(1)); 
            writeVideo(vid, F);
        end
        figure(2)
        %pcolor(x/D,z/D,u+cwave);
        pcolor(x/D,z/D,u);
        titlestr = strcat('time = ', num2str(time_vec(n)*cwave/D), ' t*c/H');
        title(titlestr);
        ylim([-1, -.85]);
        %xlim([5, 12]);
        colorbar;
        %clim([-15 15]);
        a=colorbar;
        % a.Label.String = 'u + c';
        a.Label.String = 'u';
        colormap jet;
        shading flat;
        xlabel('x/H');
        ylabel('z/H');

        drawnow;
        if(SAVEMOVIE)
            F = getframe(figure(2)); 
            writeVideo(vid2, F);
        end
        figure(3)
        pcolor(x/D,z/D,curlz);
        titlestr = strcat('time = ', num2str(time_vec(n)*cwave/D), ' t*c/H');
        title(titlestr);
        colorbar;
        ylim([-1, -.85]);
        %xlim([5, 12]);
        %clim([-15 15]);
        a=colorbar;
        a.Label.String = 'vorticity';
        colormap jet;
        shading flat;
        xlabel('x/H');
        ylabel('z/H');

        drawnow;
        if(SAVEMOVIE)
            F = getframe(figure(3)); 
            writeVideo(vid3, F);
        end
        %contourf(x/1000, z, s)
    end

    if(ENERGYPLOT)
        [r, c] = size(s);
        density = s*rho0 + rho0;
        density_vec = reshape(density, r*c, 1);
        density_sort = sort(density_vec); %should be in order of increasing salinity 
        % (increasing density)
        vol_cell = dz(1)*dx; %once dz varies, will need to edit 

        dz_background = vol_cell/(L); %flatten each density layer into 
        % a layer filling the length and width of the box, and get height 
        % of each layer. z of each will be center, so divide by 2. 

        %density should increase with depth, so define zstar vector from 0
        %to negative d 

        zstar = -dz_background/2:(-dz_background):(-D+dz_background/2);
        zstar = zstar';
        zstar = zstar+D;


        Eb(n) = sum(density_sort*g.*zstar*L*dz_background);

        z_ascend = z+D+dz(1)/2; %shift to cell centers

        Ep(n) = sum(sum(density.*g.*z_ascend*dx*dz(1)));

        Ke(n) = .5*sum(sum(density.*(u+1.71308).^2+w.^2*dx*dz(1))); 
        %Ke(n) = .5*sum(sum(density.*((u).^2+w.^2)*dx*dz(1))); 

        % figure(1)
        % pcolor(cat(2, 0*ones(length(zstar), 1), 5000*ones(length(zstar), 1)), cat(2, zstar, zstar), cat(2, density_sort, density_sort));
        % titlestr = strcat('time = ', num2str(time_vec(n)), ' seconds');
        % title(titlestr);
        % colorbar;
        % a=colorbar;
        % a.Label.String = 'density';
        % colormap jet;
        % shading flat;
        % xlabel('x (km)');
        % ylabel('z (m)');
        % 
        % drawnow;

        if(SAVEMOVIE)
            F = getframe(figure(1)); 
            writeVideo(vid, F);
        end

    end

end


if(SAVEMOVIE)
    close(vid)
    close(vid2)
    close(vid3)
end
