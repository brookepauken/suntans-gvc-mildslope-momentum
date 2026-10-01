%

% Brooke Pauken
%

close all; clear all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

exampledir = '/Users/brookepauken/suntans-gvc-mild-slope-2025/examples/quick-isw-test/mfiles/';


cd(exampledir);
ANIMATE=true;
ENERGYPLOT = true;
SAVEMOVIE=true;
SAVEUBOT=true; %make true to save figures and bottom velocities, will plot either way
TIDE=0;
initbot=0;

cwave=0;
cwave=1.75901;

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';

if ~exist('sh_val','var')
    sh_val = 2;
    nk = 200;
end
nx = 800;
dt = 0.9121;
sh_val = 2.5;

%set up datadir name and naming convention for figures
multiplier = sh_val/10;

compare = 3; %1 for dzs, 2 for xs

if compare==1
    Ns = [60, 100, 160, 200, 240, 300];
elseif compare==2
    Ns = [60, 100, 160, 200, 300, 400, 512, 600, 800, 1000];
else
    Ns = [1.75, 1, 0.6, 0.4, 0.2, 0.1];
    %Ns = [1.75, 1, 0.4, 0.1, 0.08, 0.05];
end

Nkcnt=0;

amps_v_Ns = zeros(1, length(Ns));
Ep_v_Ns = zeros(1, length(Ns));
Ke_v_Ns = zeros(1, length(Ns));
Eb_v_Ns = zeros(1, length(Ns));


for Nk = Ns
    %datadir='~/oscitests/copy_oscitest_Nz160800';
    if compare==1
        datadir=['~/oscitests/copy_oscitest_Nz', num2str(Nk), '800'];
    elseif compare==2
        datadir=['~/oscitests/copy_oscitest_Nx', num2str(Nk)];
    else
        datadir=['~/oscitests/copy_oscitest_dt', num2str(Nk)];
        %datadir=['~/oscitests/copy_oscitest_Nx400_dt', num2str(Nk)];
    end

    Nkcnt = Nkcnt+1;
    suntansfile=[datadir,'/suntans.dat'];
    
    nsteps = getvalue(suntansfile,'nsteps');
    ntout = getvalue(suntansfile,'ntout');
    dt = getvalue(suntansfile,'dt');
    nout = 1+floor(nsteps/ntout);
    
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
    D = sum(dz);
    
    
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
    
    Eb = zeros(1, nout);
    Ep = zeros(1, nout);
    Ke = zeros(1, nout);
    
    cnt = 0;
    
    ssum = zeros(1, nout);
    eta_mid=zeros(1, nout);
    max_ubot=zeros(1, nout);
    
    %middle5 = find(xv>=15000 & xv <= 35000);
    %middle5=1:5000;
    middle5=1:length(xv);
    if(initbot)
        eachside=wavel*5;
        middle5 = find(xv>=(5000-eachside-dx/2) & xv <= (5000+eachside));
        %middle5 = find(xv>=(5000-eachside) & xv <= (5000+eachside));
    end
    
    if TIDE==1
        tidal_period = (2*pi/1.45E-4); %12 hour period for M2 tide 
        num_out_tide = round(tidal_period/(ntout*dt))+1;
    
    end
    
    umid_nearbot=zeros(1, nout);
    umid_neartop=zeros(1, nout);
    
    
    %for n=1:10:(nout-10)
    for n=[1, nout]
    %for n=[1, 181]
        cnt = cnt+1;
        %fprintf('On %d of %d\n',n,nout);
        [x,z, dv, s]=plotslice_moving('s',datadir,n);   
        [~,~, ~, T]=plotslice_moving('T',datadir,n);   
        [~,~, ~, u]=plotslice_moving('u',datadir,n); 
        [~,~, ~, w]=plotslice_moving('w',datadir,n); 
        [~,~, ~, h]=plotslice_moving('h',datadir,n); 
        [~, ~, ~, q]=plotslice_moving('q', datadir,n);
    
        ubot_max(n) = max(abs(u(end, :)));
    
        if n == 1
            %D = abs(z(end, 1));
            % mean_z = multiplier*D; %height of pycnocline for calculating amplitudes
            conts = contourc(x(1, :), z(:, 1), s, 1);
            mean_z = mean(conts(2, 2:50));
            salin_thresh = conts(1, 1);
            u_base = u;
            u_old = u;
            w_base = w;
            s_base = s;
            s_old = s;
    
            h1 = multiplier*D; %0.005*D; %0.25*D;  (.1)       % Height of pycnocline
            Delta_rho = 6; %2;       % Density difference
            rho0 = 1000;         % Reference density
            g = 9.81;            % Gravity
            delta = 0.15*D; %0.1*D;       
            rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))+rho0); 
            rho_base = rho(z);
    
            z_min = min(min(z));
            z_max = max(max(z(end, :)));
    
            if(TIDE | initbot)
                bot_og = z(end, middle5);
            end
    
                [rows, cols] = size(s);
        zvec = z(:, 1)';
    
        amps_max_xs = zeros(1, cols);
    
        for ival = 1:cols
            zval = interp1(s(:, ival), zvec, salin_thresh);
            amps_max_xs(ival) = zval;
        end
    
    
    
        zmean = max(amps_max_xs); %logically this should maybe be h1 instead...
        
    
        if sh_val < 5
            zmean = max(amps_max_xs);
            %zmean = sh_val/10*D;
        else
            zmean = min(amps_max_xs);
            %zmean = sh_val/10*D;
        end
        % 
        amps_max_xs = amps_max_xs - zmean;
    
            
        end
    
            [rows, cols] = size(s);
        zvec = z(:, 1)';
    
        amps_max_xs = zeros(1, cols);
    
        for ival = 1:cols
            zval = interp1(s(:, ival), zvec, salin_thresh);
            amps_max_xs(ival) = zval;
        end
    
    
    
        zmean = max(amps_max_xs); %logically this should maybe be h1 instead...
    
    
        if sh_val < 5
            zmean = max(amps_max_xs);
            %zmean = sh_val/10*D;
        else
            zmean = min(amps_max_xs);
            %zmean = sh_val/10*D;
        end
        % 
        amps_max_xs = amps_max_xs - zmean;
    
        amps_max(n) = max(abs(amps_max_xs));
        %amps_max(n) = prctile(abs(amps_max_xs), 95);
        % 
        width(n) = (1/amps_max(n))*trapz(amps_max_xs)*dx; %characteristic width
    
    
    
        if n==1
            running_sum = u;
        else
            running_sum = running_sum + u;
        end
        
        av_u = running_sum/n;
        u_old = u;
    
        
        error(n) = sqrt(sum(sum((u-u_base).^2))/numel(u)); %try without abs
        w_error(n) = sqrt(sum(sum((w-w_base).^2))/numel(w)); %try without abs
        max_error(n) = max(max(abs(u-u_base)));
        
    
        ssum(n) = sum(sum(s));
    
        eta_mid(n)=D+z(end, length(xv)/2);
    
    
        
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
            if(length(zstar)<length(density_sort))
                zstar = [zstar, -D+dz_background/2]; %sometimes bc of rounding wont catch last point
            end
            zstar = zstar';
            zstar = zstar+D;
    
    
            Eb(n) = sum(density_sort*g.*zstar*L*dz_background);
    
            z_ascend = z+D+dz(1)/2; %shift to cell centers
    
            Ep(n) = sum(sum(density.*g.*z_ascend*dx*dz(1)));
    
            Ke(n) = .5*sum(sum(density.*(u+1.71308).^2+w.^2*dx*dz(1))); 
            %Ke(n) = .5*sum(sum(density.*((u).^2+w.^2)*dx*dz(1))); 
    
    
        end
    
        
    
        s_old = s;
        s_RMSE(n) = sqrt(sum(sum((s-s_base).^2)));
    
    end

    amps_v_Ns(Nkcnt) = amps_max(end);
    Ep_v_Ns(Nkcnt) = Ep(end);
    Ke_v_Ns(Nkcnt) = Ke(end);
    Eb_v_Ns(Nkcnt) = Eb(end);
end
%%
figure(1)
clf;
if compare==1
   total=300;
   xplot = total./Ns;
elseif compare==2
   total=10000;
   xplot = total./Ns;
else
   xplot = Ns;
   %xplot = cwave*Ns/dx;
end
%loglog(total./Ns, abs(amps_v_Ns-amps_max(1))/amps_max(1), '-o')
loglog(xplot, abs(amps_v_Ns-amps_max(1)), '-o')
hold on
%loglog(300./Ns, abs(amps_v_Ns-amps_v_Ns(end)), '-o')
title("Amplitude")
% ylabel('$|A_{5T}-A_0|/A_0$')
ylabel('$|A_{5T}-A_0|$')

hold on
Nvec = linspace(min(xplot), max(xplot), 10);
if(compare==1)
    loglog(Nvec, (.001)*(Nvec.^3), 'k--');
    xlabel('$\Delta z$')
    legend('Results', 'Slope=3', 'Location', 'best')
elseif(compare==2)
    loglog(Nvec, 1E-4*(Nvec.^2), 'k--');
    legend('Results', 'Slope=2', 'Location', 'best')
    xlabel('$\Delta x$')
else
    % load('dt800amps.mat');
    % hold on
    % plot(xplot800, yplot800, '-o')
    loglog(Nvec, 20*(Nvec.^2), 'k--');
    xlabel('$\Delta t$')  
    %xlabel('ISW Courant Number') 
    legend('Results', 'Slope=2', 'Location', 'best')

    %legend('Results Nx=400', 'Results Nx=800', 'Slope=2', 'Location', 'best')
end

%plot(linspace(min(Ns), max(Ns), 10), amps_max(1)*ones(1, 10), 'k--')
grid on


xfit = xplot(2:4);
xfit = xfit';
zfit = abs(amps_v_Ns(2:4)-amps_max(1));
zfit = zfit';
mdl = fittype('a*xfit^m','indep','xfit');
fittedmdl = fit(xfit,zfit,mdl)


% figure(2)
% semilogx(Ns, Ep_v_Ns, '-o')
% title("Potential Energy")
% grid on
% 
% figure(3)
% semilogx(Ns, Ke_v_Ns, '-o')
% title("Kinetic Energy")
% grid on
% 
figure(2)
clf;
Ea_Ns = Ep_v_Ns-Eb_v_Ns;
Ea_1 = (Ep(1)-Eb(1));
Eb_1 = Eb(1);

semilogx(xplot, (Eb_v_Ns-Eb_1)/Ea_1, '-o')
title("Available Potential Energy")
%ylabel("$(\textrm{APE}_{5T}-\textrm{APE}_0)/\textrm{APE}_0$")
ylabel("$\Delta E_b/\textrm{APE}_0$")
if(compare==1)
    xlabel('$\Delta z$')
elseif(compare==2)
    xlabel('$\Delta x$')
else
   xlabel('$\Delta t$')
   %xlabel('ISW Courant Number') 
end

grid on



% xplot800=xplot;
% yplot800 = abs(amps_v_Ns-amps_max(1));
% save('dt800amps.mat', 'xplot800', 'yplot800');