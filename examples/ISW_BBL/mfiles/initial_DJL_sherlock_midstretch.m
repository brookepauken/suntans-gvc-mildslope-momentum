%
% WRITE_INITIAL_DJL Creates a single solitary wave on the grid loaded from SUNTANS. The DJL solution is computed with
% DJLES_1.11 from Dunphy et al. (2011). Data is output for each processor into initial condition files for 
% u, s, T, and h. The u that is output is the dot product of the velocity vector with the normal on each 
% grid face. The equation of state assumes rho/rho0 = beta*s, where beta=1 in suntans.  The initial temperature field
% is a vertical column with T=1 to study the influence of mean internal wave transport
% on the scalar. Note that the minimum timestep for the internal courant number =
% 0.25 is printed but is not directly read into SUNTANS.
%

%
% Brooke Pauken
% Stanford University
%
% Heavily adapted from initialize_internal_wave.m 


if exist('eta','var') %if rerunning, clear all and set an sh_val
    clear all; 
    sh_val = 2.5;
    Atarget = 1E8;
end

%sh_val = 2.5;
%Atarget = 1E8;
%Atarget = 1.5E8;

datadir = ['../', datadir];

addpath 'DJLES-master';
addpath '/Users/brookepauken/suntans/suntans-master/mfiles';
multiplier = sh_val/10; %set fraction of depth for pycnocline

% The location of the suntans grid and files that will be written with this mfile.
suntansfile=[datadir,'/suntans.dat'];

% Write data to files. Otherwise just view the initial condition
WRITE=true;


% Output files that will be read into suntans with the processor appended to the
% end, e.g. u_init.dat.0
u_init_file = [datadir,'/u_init.dat'];
w_init_file = [datadir,'/w_init.dat'];
u_t1_file = [datadir,'/u_t1.dat'];
u_t2_file = [datadir,'/u_t2.dat'];
s_t1_file = [datadir,'/s_t1.dat'];
s_t2_file = [datadir,'/s_t2.dat'];
w_t1_file = [datadir,'/w_t1.dat'];
w_t2_file = [datadir,'/w_t2.dat'];
q_init_file = [datadir,'/q_init.dat'];
s_init_file = [datadir,'/s_init.dat'];
T_init_file = [datadir,'/T_init.dat'];
fs_init_file = [datadir,'/fs_init.dat'];
dzz_init_file = [datadir,'/dzz_init.dat'];
dzz_t1_init_file = [datadir,'/dzz_init_t1.dat'];
dzz_t2_init_file = [datadir,'/dzz_init_t2.dat'];
dzz_t3_init_file = [datadir,'/dzz_init_t3.dat'];

% Load the grid to determine the length of the domain and the cell centers
p = load([datadir,'/points.dat']);
xp = p(:,1);
L = max(xp);
cells = load([datadir,'/cells.dat']);
xv = cells(:,2);
yv = cells(:,3);

% Vertical grid
if(~exist([datadir,'/vertspace.dat']))
    error(sprintf('Can''t open %s. Try running suntans with -g first.',[datadir,'/vertspace.dat']));
end

coordsys = getvalue(suntansfile,'vertcoord');
botLayers = getvalue(suntansfile,'botLayers');
midLayers = getvalue(suntansfile,'midLayers');
cuttoffDepth = getvalue(suntansfile,'cuttoffDepth');

if(coordsys==2 || coordsys==4)
    iso =1;
else
    iso = 0;
end

if(iso) 
    dz_st = load([datadir,'/vertspace.dat']);
    dz_st = dz_st(:);
    Nk = getvalue(suntansfile,'Nkmax');
    nk_cont = Nk;
    z_st = zeros(1,Nk);
    z_st(1) = -dz_st(1)/2;
    
    %Nk = 200;
    Nk = 3000; %3000, %3800, was 6000
    dz_st = 300/Nk;
    z_st = -dz_st/2:-dz_st:(-300+dz_st/2);

    dz_st = ones(1, Nk)*dz_st;

    % higher_modes.m needs the grid to increase with increasing index, unlike the
    % suntans grid, which is the opposite
    dz_st = dz_st(end:-1:1);
    z_st = z_st(end:-1:1);

    % add bottom layers BL version w Nkbot = 45
    %mindz=0.2; %.5
    %rbot=1.0478; %1.3692

    %most sandwaves
    mindz = .2;
    Nkbot = 30;
    rbot = 1.0951;

    %higher res sandwaves
    mindz = 0.05;
    Nkbot = 40;
    rbot = 1.0963;

    %underresolved case 
    %Nkbot = 2;
    %Nkbot = 0;

    Nkbot = botLayers;
    Nkbot = 45;
    mindz = 0.2;
    rbot = 1.04878;
    if(Nkbot>0)
        dz_st_bot = zeros(1, Nkbot);
        dz_st_bot(Nkbot) = mindz;
        for k=Nkbot-1:-1:1 
            dz_st_bot(k)=rbot*dz_st_bot(k+1);
        end
        dz_st_bot=dz_st_bot(end:-1:1);

        %for constant bottom layers
        %dz_st_bot(:)=15;
        %dz_st_bot(:)=7.5;

        %Nkbot=0;
    end

else
    dz_st = load([datadir,'/vertspace.dat']);
    dz_st = dz_st(:);
    Nk = length(dz_st);
    z_st = zeros(1,Nk);
    z_st(1) = -dz_st(1)/2;

    % higher_modes.m needs the grid to increase with increasing index, unlike the
    % suntans grid, which is the opposite
    %dz_st = dz_st(end:-1:1);
    

    for k=2:Nk
        z_st(k) = z_st(k-1) - 0.5*(dz_st(k-1)+dz_st(k));
    end

    z_st = z_st(end:-1:1);
end




% Stratification is given by a tanh distribution
D = round(sum(dz_st));         % Depth
h1 = multiplier*D; %0.005*D; %0.25*D;  (.1)       % Height of pycnocline
Delta_rho = 6; %6;       % Density difference
rho0 = 1000;         % Reference density
g = 9.81;            % Gravity
delta = 0.15*D; %0.1*D;       % Pycnocline thickness

%test new strat
%delta = 25; %instead of 15
%lin_start = -150; %was -180;

lin_start = -180; %was -180;
lin_start = lin_start/300*D;

%lin_start = -310; %this should make it so there is no added strat

rho_bar = Delta_rho/2*(1-tanh((z_st+h1)/delta)); 
if(1) %temporarily turn off
    %lin_start = -180; %pick something with a slight slope still
    delta_lin = Delta_rho/2000;
    rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))); 
    rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2);
    
    %continue slope below lin_start down
    m_lin = rhoz(lin_start);
    rho_lin = m_lin*(z_st-lin_start) + rho(lin_start);
    %rho_bar(z_st<lin_start) = rho_bar(z_st<lin_start)+(lin_start-z_st(z_st<lin_start))*delta_lin;
    rho_bar(z_st<lin_start) = rho_lin(z_st<lin_start);
end

% % Density gradient and buoyancy frequency N
% drho_dz = zeros(Nk,1);
%rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2).*(z>=lin_start) + (m_lin).* (z<lin_start);
% 
drho_dz = rhoz(z_st);
% 
% % drho_dz(1) = 1.5*drho_dz(2)-0.5*drho_dz(3);
% % drho_dz(Nk) = 1.5*drho_dz(Nk-1)-0.5*drho_dz(Nk-2);
N = sqrt(-g/rho0*drho_dz);
% 
% % Compute the first mode and output the eigenfunction into phi and the
% % first-mode phase speed into c
% num_modes = 1;
% omega = 0; % Hydrostatic
% boundary_condition = 'rigid lid';
% [phi,c_linear]=higher_modes(z_st(:),N(:),omega,num_modes,boundary_condition);

%figure;
%plot(rho_bar, z_st, 'LineWidth', 1.5)
%ylabel('z (m)', 'FontSize', 20)
%xlabel('background density - $\rho_0$', 'FontSize', 20)
%xlim([0, 7]);
% lin_start = -h1*3;
% rho_bar(z_st<lin_start)
%grid on
%set(gca,'fontsize', 20) 
%pbaspect([1 1.75 1])

H = 300;
mindz = .1;

newr = Findr_Nx_minx(H, 300, mindz)

%% 

Nx_goal = length(xv);

%%%% Added in DJL code
%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Atarget = 2E8; % APE for wave (m^4/s^2) %read in as varaible
H  = D;  % domain depth (m)
Nx = Nx_goal;

%set initial Nx and Nz to be 1/4 of desired output
NX = Nx_goal/4;
if (dz_st(end)~=dz_st(1))
    NZ = (Nk*12)/4;
    too_small = 1;
    Nkmult = 12;
    cnt=0;
    if D/(Nk*Nkmult) > dz_st(end)/2
       while too_small
            Nkmult = Nkmult+4;
            if D/(Nk*Nkmult) <= dz_st(end)/2
	      too_small = 0;
              NZ = (Nk*Nkmult/4);
            end
            if cnt>10
	      too_small = 0;
            end
        end
    end
else
    NZ = Nk/4;
end
if (floor(NX) ~= NX) || (floor(Nk/4) ~= NZ)
    %if mean(dz_st)==dz_st(1)
    if dz_st(end)==dz_st(end-1)
       msg = 'Error: Number of cells in X and/or Z not divisible by 4';
       error(msg)
    end 
end

%set anon functions for rho and drho/dz
if(0)
    rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))+rho0); 
    rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2);
else
    rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))+rho0).*(z>=lin_start) + (m_lin*(z-lin_start) + rho(lin_start) + rho0).* (z<lin_start); 
    rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2).*(z>=lin_start) + (m_lin).* (z<lin_start);
end

% The velocity profile (zero for this case) (m/s) (I think we also want 0?)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;
max_iterate = 5000;
% Start with 1% of target APE, raise to target APE in 5 steps
for A=linspace(Atarget/100, Atarget, 5)
    % Find the solution of the DJL equation
    verbose = 0; %1 for full output, 2 for limited, 0 for none
    djles_refine_solution
end

% Increase the resolution, reduce epsilon, and iterate to convergence
NX=NX*2; NZ=NZ*2; epsilon=1e-6; relax = .3; max_iterate = 5000;

% Increase the resolution, and iterate to convergence
NX = NX*2; %actual grid
NZ = NZ*2; %actual grid
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
%djles_plot

djles_pressure

wavelength = djles_wavelength(eta,L)

% figure;
% pcolor(XC, ZC, pnh)
% shading flat
% colorbar



%%%% DJL code ends %%%%

%shift key outputs to edges using built in DJLES function
density_edge = djles_shift_grid(density, NX, NZ-1, 'even', 'even');
w_cent_alt = djles_shift_grid(w, NX-1, NZ, 'even', 'even');
%w_cent_alt = djles_shift_grid(wrand, NX-1, NZ, 'even', 'even');

ze_w = ze;
uwave_edge = djles_shift_grid(uwave, NX, NZ-1, 'even', 'even');
w_edge = djles_shift_grid(w, NX, NZ-1, 'heven', 'even');


xe_full = xe; %save x's of edges

%set up grid of desired suntans gridpoints and DJL gridpoints. These will differ if grid is stretched 
[XC_st, Z_st] = meshgrid(xc, z_st);
[X_edge, Z_edge] = meshgrid(xe, z_st);
[X_wst, Z_wst] = meshgrid(xc, ze);
[X_edge_DJLES, Z_edge_DJLES] = meshgrid(xe, z);

%redefine variables for clarity
ucent = uwave;
wcent = w;
uedge = uwave_edge;
wedge = w_edge;


%deal with z's being different
% Z is the output grid from DJLES, z is our true z 


if dz_st(end)~=dz_st(1)
  ucent = interp2(XC, Z, ucent, XC_st, Z_st, 'linear');
  wcent = interp2(XC, Z, wcent, XC_st, Z_st, 'linear');
  density = interp2(XC, Z, density, XC_st, Z_st, 'linear');

  uedge = interp2(X_edge_DJLES, Z_edge_DJLES, uedge, X_edge, Z_edge, 'linear');
  wedge = interp2(X_edge_DJLES, Z_edge_DJLES, wedge, X_edge, Z_edge, 'linear');
end

worse = 0;
if worse == 1 
  ucent_1 = ucent;
  wcent_1 = wcent;
  uedge_1 = uedge;
  wedge_1 = wedge;
  density_1 = density;
  XC_1 = XC;
  Z_1 = Z;
    
  %%%% run DJL code again but with larger resolution %%%%
  start_time = clock;
  max_iterate = 5000;
  clear('eta', 'var', 'ucent', 'var', 'wcent','var', 'uedge', 'var', 'wedge', 'var', 'density', 'var');
    
  NZ = (round(D/dz_st(end)))/4;
  NX = NX/4;
  % Start with 1% of target APE, raise to target APE in 5 steps
   for A=linspace(Atarget/100, Atarget, 5)
        % Find the solution of the DJL equation
        verbose = 0; %1 for full output, 2 for limited, 0 for none
        djles_refine_solution
    end
    
  % Increase the resolution, reduce epsilon, and iterate to convergence
  NX=NX*2; NZ=NZ*2; epsilon=1e-6; relax = .3; max_iterate = 5000;
    
% Increase the resolution, and iterate to convergence
    NX = NX*2; %actual grid
    NZ = NZ*2; %actual grid
    djles_refine_solution
    end_time=clock;
    fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));
    
    % Compute and plot the diagnostics
    djles_diagnostics
    djles_plot
    
    %%%% 2nd DJL code ends %%%%
    
    %shift key outputs to edges using built in DJLES function
    density_edge_2 = djles_shift_grid(density, NX, NZ-1, 'even', 'even');
    uwave_edge_2 = djles_shift_grid(uwave, NX, NZ-1, 'even', 'even');
    w_edge_2 = djles_shift_grid(w, NX, NZ-1, 'even', 'even');
    
    %xe_full = xe; %save x's of edges
    
    %set up grid of desired suntans gridpoints and DJL gridpoints. These will differ if grid is stretched 
    [X_edge_DJLES_2, Z_edge_DJLES_2] = meshgrid(xe, z);
    
    %redefine variables for clarity
    ucent_2 = uwave;
    wcent_2 = w;
    uedge_2 = uwave_edge;
    wedge_2 = w_edge;


    %interpolate finer solution onto grid 
    ucent_1 = interp2(XC_1, Z_1, ucent_1, XC_st, Z_st, 'linear');
    wcent_1 = interp2(XC_1, Z_1, wcent_1, XC_st, Z_st, 'linear');
    density_1 = interp2(XC_1, Z_1, density_1, XC_st, Z_st, 'linear');

    uedge_1 = interp2(X_edge_DJLES, Z_edge_DJLES, uedge_1, X_edge, Z_edge, 'linear');
    wedge_1 = interp2(X_edge_DJLES, Z_edge_DJLES, wedge_1, X_edge, Z_edge, 'linear');

    %find cutoff point
    z_coarse_inds = find(dz_st==dz_st(end));
    z_cutoff = z_st(z_coarse_inds(1));
    z_coarse_inds_2 = find(z>=z_cutoff);
    if(length(z_coarse_inds)>length(z_coarse_inds_2))
        z_coarse_inds_2 = [z_coarse_inds_2-1, z_coarse_inds_2(end)];
    end

    ucent_1(z_coarse_inds(1), :) = .5*(ucent_2(z_coarse_inds_2(1), :) + ucent_1(z_coarse_inds(1), :));
    wcent_1(z_coarse_inds(1), :) = .5*(wcent_2(z_coarse_inds_2(1), :) + wcent_1(z_coarse_inds(1), :));
    density_1(z_coarse_inds(1), :) = .5*(density(z_coarse_inds_2(1), :) + density_1(z_coarse_inds(1), :));

    uedge_1(z_coarse_inds(1), :) = .5*(uedge_2(z_coarse_inds_2(1), :)+uedge_1(z_coarse_inds(1), :));
    wedge_1(z_coarse_inds(1), :) = .5*(wedge_2(z_coarse_inds_2(1), :)+wedge_1(z_coarse_inds(1), :));
    
    %replace top section with coarser grid solution 
    ucent_1(z_coarse_inds(2:end), :) = ucent_2(z_coarse_inds_2(2:end), :);
    wcent_1(z_coarse_inds(2:end), :) = wcent_2(z_coarse_inds_2(2:end), :);
    density_1(z_coarse_inds(2:end), :) = density(z_coarse_inds_2(2:end), :);

    uedge_1(z_coarse_inds(2:end), :) = uedge_2(z_coarse_inds_2(2:end), :);
    wedge_1(z_coarse_inds(2:end), :) = wedge_2(z_coarse_inds_2(2:end), :);

    %set the first coarse value to average of the two grids 


    ucent=ucent_1;
    wcent=wcent_1;
    density=density_1;
    uedge=uedge_1;
    wedge=wedge_1;

    clear('ucent_1', 'var', 'ucent_2', 'var', 'wcent_1', 'var', 'wcent_2', 'var', 'density_1', 'var', 'density_2', 'var', 'uedge_1', 'var', 'uedge_2', 'var', 'wedge_1', 'var', 'wedge_2', 'var');


end

figure;
pcolor(XC_st, Z_st, density)
shading flat
colorbar


% c_linear / c

%%
% figure;
% pcolor(XC/1000, Z, ucent)
% shading flat
% a = colorbar
% a.Label.Interpreter = 'latex';
% a.Label.FontSize=12;
% a.Label.String = 'u $(m/s)$';
% colormap jet;
% shading flat;
% xlabel('x (km)', 'FontSize', 12);
% ylabel('z (m)', 'FontSize', 12);
% 
% figure;
% pcolor(X_edge_DJLES/1000, Z_edge_DJLES, uedge)
% shading flat
% a = colorbar
% a.Label.Interpreter = 'latex';
% a.Label.FontSize=12;
% a.Label.String = 'uedge $(m/s)$';
% colormap jet;
% shading flat;
% xlabel('x (km)', 'FontSize', 12);
% ylabel('z (m)', 'FontSize', 12);
% 
% 
% c_linear / c



xc_st = xc;
xe_st = xe;


if(iso) %%if isopycnal, need to write layer heights, etc
    isoloopstart = tic;
    parpool('local', 8) %start parallel loop with 8 processors
    %parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))

    ncent = length(xc_st);
    x_edge_list = X_edge(1, :);
    nedge = length(xe_st);

    dx_st = L/Nx;

    % 
    % layers_add_bottom = Nkbot; %4
    % layers_add_mid = 10; %was 15, testing double 10 (10), 20 for normal resolved bedforms!!
    % layers_add_top = 0; %2 (4)

    %randomness max and min
    a = 0;
    b = 0;

    
    %normal iso
    layers_add_bottom = Nkbot; %4
    layers_add_mid = midLayers; %This is not always -2 but it is a good starting point.
    layers_add_mid_iso = 0;
    layers_add_top = 0; %4

    mainloopmax = (nk_cont-layers_add_bottom);

    u_layer_av_cent = zeros(nk_cont-(layers_add_bottom+layers_add_mid), length(xc_st));
    w_layer_av_cent = zeros(nk_cont-(layers_add_bottom+layers_add_mid), length(xc_st));
    density_layer_av = zeros(nk_cont-(layers_add_bottom+layers_add_mid), length(xc_st));
    layer_thickness = zeros(nk_cont-(layers_add_bottom+layers_add_mid), length(xc_st));

    top_dens = density(1, 1);
    bot_dens = density(end, end);
    
    
    
    x_cont_lower = x_edge_list;
    if(layers_add_bottom>0)
        relative_bot = 270; %270 for most sand wave runs. 280 for 90 runs
        relative_bot = H-sum(dz_st_bot);
        buffer = 3;
    else
        relative_bot = H;
        buffer = 0; %1 maybe
    end

    z_cont_lower = 0*x_cont_lower-relative_bot;
    
    
    x_cont_top = x_edge_list;
    z_cont_top = 0*x_cont_top;
    
    XE_list = X_edge(:);
    ZE_list = Z_edge(:);
    
    XC_list = XC(:);
    ZC_list = ZC(:);

    figure;

    % First do contour with Ntop+1 layers (we wont use the bottom one,
    % we'll replace it with the middle and bottom). Write these layers as
    % usual, but probably change density definition 
    
    % want one extra layer 
    % no we dont?
    contour(XC(1, :), Z(:, 1), density, nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top), 'k') %removed +1
    %conts = contourc(XC(1, :), Z(:, 1), density, nk_cont-1-layers_add_bottom);
    conts = contourc(XC(1, :), Z(:, 1), density, nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top));
    salvals = zeros(1, nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top));
    
    if(layers_add_top>0)                                                                                                                                             
        identifier_ind = 1;                                                                                                                                          
                                                                                                                                                                     
        %remove isopycnals some depth                                                                                                                               
                                                                                                                                                                     
         cutoff = -15; %-20                                                                                                                                              
         remove_layers_top = 0;    
         first_zs = zeros(1, (nk_cont-(layers_add_bottom+layers_add_mid+layers_add_top)));                                                                          
                                                                                                                                                                     
        for layercnt=1:(nk_cont-(layers_add_bottom+layers_add_mid+layers_add_top))                                                                                
           last_ind = identifier_ind+conts(2, identifier_ind);                                                                                                      
           salvals(layercnt) = conts(1, identifier_ind);                                                                                                            
            first_zs(layercnt) = conts(2, identifier_ind+1);                                                                                                                   
                                                                                                                                                                     
                                                                                                                                                                     
           if(first_zs(layercnt) > cutoff)                                                                                                                                 
             remove_layers_top = remove_layers_top+1;                                                                                                                
           end                                                                                                                                                       
           identifier_ind = last_ind+1;                                                                                                                             
    end
                                                                                                                                                                                                                                                                                                                           
        if(remove_layers_top>0)
        highest_layer_height = first_zs(remove_layers_top+1)-first_zs(remove_layers_top+2);
    else
        highest_layer_height = first_zs(1)-first_zs(2);
    end

        layers_add_top = layers_add_top+remove_layers_top;

    zs_tops = linspace(first_zs(remove_layers_top+1), 0, layers_add_top+2);
    zs_tops = zs_tops(2:end-1);
    s_tops = interp1(z_st, density(:, 1),  zs_tops);

    %layers_add_top = layers_add_top + remove_layers_top;                                                                                                         
    salvals(1:remove_layers_top) = [];                                                                                                                           
                                                                                                                                                                    
    %add_salvals = linspace(min(min(density)), min(salvals), layers_add_top+2);                                                                                    
    salvals = [s_tops, salvals];                                                                                                                   
    conts = contourc(XC(1, :), Z(:, 1), density, salvals);                                                                                                              

    figure;
    contour(XC(1, :), Z(:, 1), density, salvals((layers_add_top+1):end), 'k')
    xlabel('x (m)')
    ylabel('z (m)')
    hold on
    contour(XC(1, :), Z(:, 1), density, salvals(1:layers_add_top), 'b')
    else
        remove_layers_top = 0;
    end  
                             

    if(layers_add_mid_iso==1)                                                                                                                                             
      identifier_ind = 1;                                                                                                                                          
                                                                                                                                                                     
        %remove isopycnals some depth                                                                                                                               
                                                                                                                                                                     
	cutoff = -150; %-20                                                                                                                                              
	remove_layers_bot = 0;    
first_zs = zeros(1, (nk_cont-(layers_add_bottom+layers_add_mid+layers_add_top))-1);
                                                                                                                                                                     
for layercnt=1:(nk_cont-(layers_add_bottom+layers_add_mid+layers_add_top)-1)                                                                                
	       last_ind = identifier_ind+conts(2, identifier_ind);                                                                                                      
salvals(layercnt) = conts(1, identifier_ind);                                                                                                            
first_zs(layercnt) = conts(2, identifier_ind+1);                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
identifier_ind = last_ind+1;                                                                                                                             
        end                                                                                                                                                              

        lowest_layer_height = first_zs(end)-first_zs(end-1);

    
zs_bots = linspace(first_zs(end), -H, layers_add_mid+2);
zs_bots= zs_bots(2:end-1);
s_bots = interp1(z_st, density(:, 1),  zs_bots);
    
%layers_add_top = layers_add_top + remove_layers_top;  
%salvals(1:remove_layers_top) = [];                                                                                                                           
                                                                                                                                                                     
%add_salvals = linspace(min(min(density)), min(salvals), layers_add_top+2);                                                                                    
salvals = [salvals, s_bots];                                                                                                                   
conts = contourc(XC(1, :), Z(:, 1), density, salvals);  
    
figure;
contour(XC(1, :), Z(:, 1), density, salvals((layers_add_top+1):end), 'k')
xlabel('x (m)')
ylabel('z (m)')
        hold on
contour(XC(1, :), Z(:, 1), density, salvals(1:layers_add_top), 'b')

layers_add_mid = 0;
 else
   remove_layers_mid = 0;
    end  

                                                                                                                                                                     
    if(layers_add_mid>0)
    first_z = zeros(1, nk_cont-1-(layers_add_bottom+layers_add_mid));                                                                                                
    identifier_ind = 1; 
        

    %remove isopycnals some depth
    %cutoff = -159; %150 for "new" strat
        cutoff = -150; %for basic runs
    remove_layers = 0;
        maxlayer = 40; %20 (was 40!!!!)


    for layercnt=1:(length(salvals))
        % if(layercnt==length(salvals))
        %     last_ind = length(conts);
        % else
        %     last_ind = identifier_ind+conts(2, identifier_ind);
        % end

       last_ind = identifier_ind+conts(2, identifier_ind);
           

        salvals(layercnt) = conts(1, identifier_ind);
        first_z(layercnt) = conts(2, identifier_ind+1);
        if(first_z(layercnt)<cutoff || first_z(layercnt)<cutoff+maxlayer && abs(first_z(layercnt-1) - first_z(layercnt))> maxlayer)
            remove_layers = remove_layers+1;
                for extralayer=(layercnt+1):(nk_cont-1-(layers_add_bottom+layers_add_mid))
                remove_layers = remove_layers+1;
            end
            break;
        end
        old_ind = identifier_ind;
        identifier_ind = last_ind+1;
        
    end

        if(layers_add_bottom==0)
        layers_add_mid = layers_add_mid + remove_layers; %removed +1
        else
            layers_add_mid = layers_add_mid + remove_layers +1; %removed +1
        end


    z_bot_top = conts(2, (old_ind+1):(old_ind+conts(2, old_ind))); %these are the z values for the bottom of the top section
    x_bot_top = conts(1, (old_ind+1):(old_ind+conts(2, old_ind)));


    %z_bot_top_cent = interp1(x_bot_top, z_bot_top, xc); %now these are z values but at x centers
    z_bot_top_edge = interp1(x_bot_top, z_bot_top, xe); %now these are z values but at x edges

    z_add = zeros(layers_add_mid+1, length(xe));

    for xlayer = 1:length(xe)
        z_add(:, xlayer) = linspace(-relative_bot, z_bot_top_edge(xlayer), layers_add_mid+1);
    end

    z_add(:, 1) = z_add(:, 2);
    z_add(:, end) = z_add(:, end-1);
    z_add(end, :) = [];

        % z_add(:, 1) = z_add(:, 2);
        % z_add(:, end) = z_add(:, end-1);
        % z_add(end, :) = [];
        % z_add = z_add(2:end, :);

        dz_mid = zeros(layers_add_mid+1, length(xe));
    
           STRETCH_MID=1;
           if(STRETCH_MID)
                %get stretching factor from top bottom layer to top iso at
                %edge
                r = Findr_Nx_minx(z_bot_top_edge(2)--relative_bot+max(dz_st_bot), layers_add_mid+1, max(dz_st_bot));

                z_bot_top_edge(1) = z_bot_top_edge(2);
                z_bot_top_edge(end) = z_bot_top_edge(end-1);


                dz_mid(layers_add_mid+1, 1) = max(dz_st_bot);
                %dz(Nkmax) = H*(r-1)/(r^Nkmax-1);
                % dz(Nkmax) = mindz;
                for k=layers_add_mid:-1:1 
                    dz_mid(k, 1)=r*dz_mid(k+1, 1);
                end

                %now we have the dzs for the middle, we need to find dsig 
                dsig_mid = dz_mid(1:layers_add_mid, 1)./(z_bot_top_edge(2)--relative_bot);

                %now apply to the remaining distances to get the rest of
                %the dzs

                dz_mid(layers_add_mid+1, :) = max(dz_st_bot);
                dz_mid(1:layers_add_mid, :) = dsig_mid.*(z_bot_top_edge--relative_bot);

                x_mid = ones(layers_add_mid, 1)*xe;

                plot(x_mid', cumsum(dz_mid((end-1):-1:1, :))'-relative_bot, 'b');
                hold on 
                plot(xe, z_bot_top_edge, 'k')
                hold on
                plot(xe, ones(1, length(xe))*(-relative_bot), 'k');


                z_add(2:end, :) = -relative_bot + cumsum(dz_mid((end-1):-1:2, :));
                z_add(1, :) = -relative_bot;
                


            end

        if(layers_add_bottom==0)
    
	    % z_add = zeros((layers_add_mid+1)*2, length(xe)); %+1?
            % 
	    % for xlayer = 1:length(xe)
	    %     z_add(:, xlayer) = linspace(-relative_bot+buffer, z_bot_top_edge(xlayer), (layers_add_mid+1)*2);
            % end
            % 
            % z_add(:, 1) = z_add(:, 2);
            % z_add(:, end) = z_add(:, end-1);
            % 
            % z_add = z_add(2:2:end-1, :);

            z_add = zeros(layers_add_mid+2, length(xe));
    
            for xlayer = 1:length(xe)
	       z_add(:, xlayer) = linspace(-relative_bot, z_bot_top_edge(xlayer), layers_add_mid+2);
            end
        
            z_add(:, 1) = z_add(:, 2);
            z_add(:, end) = z_add(:, end-1);
            z_add(end, :) = [];
            z_add = z_add(2:end, :);
    
            s_bots = interp1(z_st, density(:, 1),  z_add(:, 1));
    
            s_bots = mean(s_bots, 2);
            s_bots = s_bots(end:-1:1);
        
            %layers_add_top = layers_add_top + remove_layers_top;  
            salvals((end-(remove_layers)):end) = [];                                                                                                                           
                                                                                                                                                                         
           %add_salvals = linspace(min(min(density)), min(salvals), layers_add_top+2);                                                                                    
            salvals = [salvals, s_bots'];                                                                                                                   
            conts = contourc(XC(1, :), Z(:, 1), density, salvals);  
        
            figure;
            contour(XC(1, :), Z(:, 1), density, salvals(1:end), 'g')
            hold on
            contour(XC(1, :), Z(:, 1), density, salvals((layers_add_top+1):end-length(s_bots)), 'k')
            xlabel('x (m)')
            ylabel('z (m)')
            hold on
            contour(XC(1, :), Z(:, 1), density, salvals(1:layers_add_top), 'b')
            hold on
    
            if(layers_add_bottom==0)
            layers_add_mid = 0;
            end
        end
        
    else
         remove_layers = 0;
         z_add = zeros(nk_cont, Nx);
    end
    % think we can skip adding more bottom density layers, but we'll see 

    % z of top of last layer is the top of the middle. For each x location
    % divide up remaining space evenly for z's. Start a z vec and we'll add
    % that to bottom section?

    % 
    % contour(XC(1, :), Z(z_st>=-relative_bot, 1), density(z_st>=-relative_bot, :), nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top), 'k')
    % %conts = contourc(XC(1, :), Z(:, 1), density, nk_cont-1-layers_add_bottom);
    % conts = contourc(XC(1, :), Z(z_st>=-relative_bot, 1), density(z_st>=-relative_bot, :), nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top));
    % salvals = zeros(1, nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top));
    % first_z = zeros(1, nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top));
    % identifier_ind = 1;
    % for layercnt=1:(nk_cont-1-(layers_add_bottom+layers_add_mid+layers_add_top))
    %     last_ind = identifier_ind+conts(2, identifier_ind);
    %     salvals(layercnt) = conts(1, identifier_ind);
    %     first_z(layercnt) = conts(2, identifier_ind+1);
    %     identifier_ind = last_ind+1;
    % end

    %first loop over and save identifier inds 
    start_inds = ones(1, (nk_cont-layers_add_bottom-layers_add_mid)+1);

    if layers_add_bottom==0
        sub1 = 1;
    else
        sub1 = 0;
    end

    for layercnt = 1:(nk_cont-layers_add_bottom-layers_add_mid-sub1) %removed minus 1 from limit then added it back
        
        start_inds(layercnt+1) = start_inds(layercnt)+conts(2, start_inds(layercnt))+1;

        % hold on
        % layercnt
        % h1 = plot(conts(1, (start_inds(layercnt)+1):start_inds(layercnt+1)-1), conts(2, (start_inds(layercnt)+1):start_inds(layercnt+1)-1), 'r');
        % pause;
        % delete(h1);

        

        % if(start_inds(layercnt)>length(conts))
        %     start_inds(layercnt+1) = length(conts);
        % else
        %     start_inds(layercnt+1) = start_inds(layercnt)+conts(2, start_inds(layercnt))+1;
        % end

    end
    if(remove_layers==0 || layers_add_bottom==0)
        if(start_inds(end-1)>length(conts))
            start_inds(end-1) = floor(length(conts)-1);
            start_inds(2:end) = start_inds(1:end-1);
        else
            start_inds(end) = length(conts)-1;
        end
        %start_inds(end) = length(conts)-1;
        
    else
        layercnt=layercnt; %removed +1
        start_inds(layercnt+1) = start_inds(layercnt)+conts(2, start_inds(layercnt))+1;
    end

    start_inds = start_inds(end:-1:1);
    start_inds = cat(2, ones(1, layers_add_mid-1), start_inds); %added -1


    identifier_ind = 1;
    s_layer = 0;
    s_cont_top = rho0;
        mainloopvals = 1:(nk_cont-layers_add_bottom);
    flipmainloopvals = (nk_cont-layers_add_bottom)-1:1;
    parfor layercnt=1:(nk_cont-layers_add_bottom)
    %for layercnt=1:(nk_cont-layers_add_bottom)

        x_cont_top = x_edge_list;
        z_cont_top = 0*x_cont_top;
        s_layer = 0;
        if((layercnt ~= 1 && layercnt>layers_add_mid) || layers_add_mid==0)

            %last_ind = identifier_ind+conts(2, identifier_ind);
            if(layercnt>1)
                identifier_ind = start_inds(layercnt)+1; %+1
                last_ind = start_inds(layercnt-1)-1; %0
                s_cont = conts(1, identifier_ind-1);
                x_cont = conts(1, identifier_ind:last_ind);
                z_cont = conts(2, identifier_ind:last_ind);
            end


            % %z_cont = conts(2, identifier_ind+1:last_ind);
            % 
            % %%add random variation to layer height
            % [rows, cols] = size(z_cont);
            % rands = a + (b-a).*rand(rows,cols);
            % z_cont = z_cont+rands;
            

            if(layercnt==(nk_cont-layers_add_bottom))
                %s_cont_top = min(min(density));
                s_cont_top = rho(0);
                x_cont_top = x_edge_list;
                z_cont_top = 0*x_cont_top;
            elseif(layercnt==1)
                s_cont_top = conts(1, start_inds(layercnt+1));
                x_cont_top = conts(1, start_inds(layercnt+1)+1:end);
                z_cont_top = conts(2, start_inds(layercnt+1)+1:end);


                x_cont_top = cat(2, 0, x_cont_top);
                x_cont_top = cat(2, x_cont_top, L);
                
                z_cont_top = cat(2, z_cont_top(1), z_cont_top);
                z_cont_top = cat(2, z_cont_top, z_cont_top(end)); 

                %s_cont = max(max(density));
                s_cont = rho(-relative_bot);
                x_cont = xe;
                x_cont(1) = [];
                x_cont(end) = [];


                z_cont = 0*x_cont - D;
                z_cont = 0*x_cont - relative_bot;
                % z_cont = cat(2, z_cont(1), z_cont);
                % z_cont = cat(2, z_cont, z_cont(end));
            else
                identifier_ind = start_inds(layercnt+1); %+2
                last_ind = start_inds(layercnt)-1; %+1
                    

                s_cont_top = conts(1, start_inds(layercnt+1)); %+2
                x_cont_top = conts(1, identifier_ind+1:last_ind);
                z_cont_top = conts(2, identifier_ind+1:last_ind);

                x_cont_top = cat(2, 0, x_cont_top);
                x_cont_top = cat(2, x_cont_top, L);
                z_cont_top = cat(2, z_cont_top(1), z_cont_top);
                z_cont_top = cat(2, z_cont_top, z_cont_top(end));
            end

            s_layer = .5*(s_cont_top + s_cont);
    
    
            x_cont = cat(2, 0, x_cont);
            x_cont = cat(2, x_cont, L);
            
            z_cont = cat(2, z_cont(1), z_cont);
            z_cont = cat(2, z_cont, z_cont(end));

        elseif(layercnt==layers_add_mid)
            x_cont = xe;
            z_cont = z_add(layercnt, :);

            identifier_ind = start_inds(layercnt+1);
            last_ind = start_inds(layercnt)-1;   


            x_cont_top = conts(1, identifier_ind+1:last_ind);
            z_cont_top = conts(2, identifier_ind+1:last_ind);

            x_cont_top = cat(2, 0, x_cont_top);
            x_cont_top = cat(2, x_cont_top, L);
            z_cont_top = cat(2, z_cont_top(1), z_cont_top);
            z_cont_top = cat(2, z_cont_top, z_cont_top(end));
            
        elseif(layercnt ~= 1 && layercnt<layers_add_mid)
            x_cont = xe;
            x_cont_top = xe;
            z_cont = z_add(layercnt, :);
            z_cont_top = z_add(layercnt+1, :);


        else %bottom
            if(layers_add_mid>0)
                x_cont_top = xe;
                z_cont_top = z_add(layercnt+1, :);

            else
                identifier_ind = start_inds(layercnt+2);
                last_ind = start_inds(layercnt+1)-1;

                x_cont_top = conts(1, identifier_ind+1:last_ind);
                z_cont_top = conts(2, identifier_ind+1:last_ind);


                x_cont_top = cat(2, 0, x_cont_top);
                x_cont_top = cat(2, x_cont_top, L);
                z_cont_top = cat(2, z_cont_top(1), z_cont_top);
                z_cont_top = cat(2, z_cont_top, z_cont_top(end));
            end

            x_cont = x_cont_lower;
            z_cont = z_cont_lower;

    
        end

        xp = cat(2, x_cont_top(end:-1:1), x_cont);
        zp = cat(2, z_cont_top(end:-1:1), z_cont);

        xp = cat(2, xp, xp(1));
        zp = cat(2, zp, zp(1));

        in = inpolygon(XE_list, ZE_list, xp', zp');
        in_center = inpolygon(XC_list, ZC_list, xp', zp');
        x_in = XE_list(in);
        z_in = ZE_list(in);

        x_in_c = XC_list(in_center);
        z_in_c = ZC_list(in_center);
        in_indicies = 1:length(XE_list);
        in_indicies = in_indicies(in);

        in_indicies_cent = 1:length(XC_list);
        in_indicies_cent = in_indicies_cent(in_center);


        % %plot to visualize if desired
        % if(layercnt<15)
        %     hold on;
        %     layerhigh = plot(xp, zp, 'r');
        %     pause;
        %     delete(layerhigh);
        % end

        % hold on;
        % layerhigh = plot(xp, zp, 'r');
        % pointshigh = scatter(x_in_c, z_in_c, 'b');
        % pause;
        % delete(layerhigh);
        % delete(pointshigh);

        %xcnt=0; 
        u_in = ucent(in_center);
        w_in = wcent(in_center);
        
        for xcnt = 1:ncent
            % xcnt = xcnt+1;
            %xlayer = xc(xcnt);
            xlayer = xc_st(xcnt);


            u_layer_av_cent(layercnt, xcnt) = mean(u_in(abs(x_in_c-xlayer)<dx_st/2));
            w_layer_av_cent(layercnt, xcnt) = mean(w_in(abs(x_in_c-xlayer)<dx_st/2));
            
            if(layercnt>layers_add_mid)
                density_layer_av(layercnt, xcnt) = s_layer;
            else
                density_in = density(in_center);
                density_layer_av(layercnt, xcnt) = mean(density_in(abs(x_in_c-xlayer)<dx_st/2));
            end
            layer_thickness(layercnt, xcnt) = interp1(x_cont_top, z_cont_top, xlayer) - interp1(x_cont, z_cont, xlayer);
            %u_layer_av_cent(layercnt, xcnt) = 1/layer_thickness(layercnt, xcnt).*sum(u_in.*(z_in_c-z_in_c(end)));
            q_layer_av(layercnt, xcnt) = mean(pnhz(abs(x_in_c-xlayer)<dx_st/2));
        end

        u_in = uedge(in);
        w_in = wedge(in);

        for xcnt = 1:nedge
            xlayer = x_edge_list(xcnt);
            xlayer = xe_st(xcnt);
            indsx = find(X_edge==xlayer);
            layer_thickness_edge(layercnt, xcnt)=interp1(x_cont_top, z_cont_top, xlayer) - interp1(x_cont, z_cont, xlayer);


            %u_layer_av(layercnt, xcnt) = 1/layer_thickness_edge.*sum(u_in.*(z_in-z_in(end)));
            u_layer_av(layercnt, xcnt) = mean(u_in(abs(x_in-xlayer)<dx_st/2));
            w_layer_av(layercnt, xcnt) = mean(w_in(abs(x_in-xlayer)<dx_st/2));
        end

        z_cont = z_cont_top;

    end

    %  add top layers if desired (put before loop above)
    % add_salvals = linspace(min(min(density)), min(salvals), layers_add_top+2);
    % salvals = [add_salvals(2:end-1), salvals];
    % conts = contourc(XC(1, :), Z(z_st>=-relative_bot, 1), density(z_st>=-relative_bot, :), salvals);
    % figure;
    % contour(XC(1, :), Z(z_st>=-relative_bot, 1), density(z_st>=-relative_bot, :), salvals, 'b')

    %do reamining layers (non density)
    if(layers_add_bottom>0)
       % z_remain = linspace(-relative_bot, -H, layers_add_bottom*2+1);
       % bot_layer_thicness = ones(layers_add_bottom, length(xc_st))*(z_remain(1)-z_remain(3));
       % bot_layer_thicness_edge = ones(layers_add_bottom, length(xe_st))*(z_remain(1)-z_remain(3));
       %for stretching 
        z_remain = -(300-cumsum(dz_st_bot));
        z_remain = [-300, z_remain];

        %make sure this is right
        z_remain(end)=-relative_bot;
        dz_st_bot(end) = z_remain(end)-z_remain(end-1);

        %now put at centers
        z_remain = .5*(z_remain(2:end)+z_remain(1:end-1));

        %for stretching 
        z_remain_w = -(300-cumsum(dz_st_bot));
        z_remain_w = [-300, z_remain_w];

        %make sure this is right
        z_remain_w(end)=-relative_bot;

        bot_layer_thicness = ones(length(xc), 1)*dz_st_bot;
        bot_layer_thicness = bot_layer_thicness';
        bot_layer_thicness_edge = ones(length(xe), 1)*dz_st_bot;
        bot_layer_thicness_edge = bot_layer_thicness_edge';
        layer_thickness = cat(1, bot_layer_thicness, layer_thickness);
        layer_thickness_edge = cat(1, bot_layer_thicness_edge, layer_thickness_edge);

        [Xe_bot, Ze_bot] = meshgrid(xe_st, z_remain);
        [Xc_bot, Zc_bot] = meshgrid(xc_st, z_remain);
        [Xc_w, Zc_w] = meshgrid(xc_st, z_remain_w);
        XC_wst = cat(1, XC_st(1, :),  XC_st);


        ucent_bot = interp2(XC, ZC, ucent, Xc_bot, Zc_bot, 'linear');
        wcent_bot = interp2(XC_wst, Z_wst, w_cent_alt, Xc_w, Zc_w, 'linear');
        density_bot = interp2(XC, ZC, density, Xc_bot, Zc_bot, 'linear');
        q_layer_bot = interp2(XC, ZC, pnhz, Xc_bot, Zc_bot, 'linear');

        uedge_bot = interp2(X_edge_DJLES, Z_edge_DJLES, uedge, Xe_bot, Ze_bot, 'linear');
        wedge_bot = interp2(X_edge_DJLES, Z_edge_DJLES, wedge, Xe_bot, Ze_bot, 'linear');


        u_layer_av = cat(1, uedge_bot(:, :), u_layer_av);
        u_layer_av_cent = cat(1, ucent_bot(:, :), u_layer_av_cent);
        w_layer_av_cent = cat(1, wcent_bot(:, :), w_layer_av_cent);
        density_layer_av = cat(1, density_bot(:, :), density_layer_av);
        q_layer_av = cat(1, q_layer_bot(:, :), q_layer_av);
    end

    % force depth averaged velocity to be 0
    % uedge = u_layer_av-(1/H.*sum(u_layer_av.*layer_thickness_edge, 1));
    % ucent = u_layer_av_cent-1/H.*sum(u_layer_av_cent.*layer_thickness, 1);

    %try forcing such that w has no depth average
    %wcent = w_layer_av_cent-(1/H.*sum(w_layer_av_cent.*layer_thickness, 1));
    %check
    %max(max(1/H.*sum(wcent.*layer_thickness, 1)))

    %origianl
    wcent = w_layer_av_cent;

    uedge = u_layer_av;
    ucent = u_layer_av_cent;
    density = density_layer_av;
    [X_edge_new, ~] = meshgrid(x_edge_list, 1:nk_cont);
    Nk = nk_cont;

    %[w_zedge, w_zcent, diffs] = w_continuity(wcent, layer_thickness, uedge, Nk, NX);
    %wcent = w_zcent; %replace with continuity version
    % figure;
    % pcolor(XC_st, cumsum(layer_thickness, 1), w_zcent);
    % colorbar
    % shading flat
    % colormap jet;
    % 
    % figure;
    % pcolor(XC_st, cumsum(layer_thickness, 1), wcent);
    % colorbar
    % shading flat
    % colormap jet;
    % 
    % figure;
    % pcolor(XC_st, cumsum(layer_thickness, 1), abs(diffs));
    % colorbar
    % shading flat
    % colormap jet;

    delete(gcp);

    Tprint = toc(isoloopstart);
    fprintf('Total parallel iso loop executed in %8.2f seconds \n', Tprint);

end



dt = getvalue(suntansfile,'dt');
dtuse = dt;

Nsteps = L/c * 3 / dt 

Courantnum = dt*c/dx

dtprint = (.01)*(dx)/c

dtprint_ualt = (.01)*(dx)/max(max(ucent))



Nsteps_20 = L/c * 20 / min(dtprint, dtprint_ualt)

if(iso)
    [XC_st, Z_st] = meshgrid(xc_st, linspace(D/(2*nk_cont), D-D/(2*nk_cont), nk_cont));
    %Z_st = cumsum(layer_thickness, 1);
    [X_edge, Z_edge] = meshgrid(xe_st, 1:nk_cont);
    Z_newgrid = cumsum(layer_thickness, 1);
    Z_newgrid_edge = cumsum(layer_thickness_edge, 1);
    % Z_newgridtemp = Z_newgrid;
    % Z_newgridtemp_edge = Z_newgrid_edge;
    % 
    % Z_newgridtemp(1, :) = layer_thickness(1, :)./2;
    % Z_newgridtemp_edge(1, :) = layer_thickness_edge(1, :)./2;
    % for k = 2:nk_cont
    %     Z_newgridtemp(k, :) = Z_newgrid(k-1, :) + layer_thickness(k, :)./2;
    %     Z_newgridtemp_edge(k, :) = Z_newgrid_edge(k-1, :) + layer_thickness_edge(k, :)./2;
    % end
    
else
    Z_newgrid = Z_st;
    Z_newgrid_edge = Z_edge;   
end

Z_newgrid = Z_st;
Z_newgrid_edge = Z_edge;  

% Z_newgrid = Z_newgridtemp;
% Z_newgrid_edge = Z_newgridtemp_edge;  
%% new section to get w on top and bottom faces

 if(iso)
   Z_w_st = cumsum(layer_thickness, 1);                                                                                                                             
   Z_w_st = cat(1, zeros(1, length(xc_st)), Z_w_st);                                                                                                                
   Z_w_st = Z_w_st + -300;
   Z_other_st = cumsum(layer_thickness, 1) - layer_thickness./2 - 300;  
   X_shift_w = cat(1, XC_st(1, :),  XC_st);
   [~, Z_shift_w] = meshgrid(xc_st, linspace(0, 300, nk_cont+1));
else
   %[~, Z_shift_w] = meshgrid(xc_st, linspace(0, 300, Nk+1));                                                                                                  
   %X_shift_w = cat(1, XC_st(1, :),  XC_st);                                                                                                                        
   Z_shift_w = Z_st;
   X_shift_w = XC_st;
end

%%
%shift everything back one and two time steps to write t0-dt and t0-2dt
ucent_t1 = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, ucent, Z_newgrid);
ucent_t2 = timeshift_soln(c, dtuse, L, 2, XC_st, Z_st, ucent, Z_newgrid);
uedge_t1 = timeshift_soln(c, dtuse, L, 1, X_edge, Z_edge, uedge, Z_newgrid_edge);
uedge_t2 = timeshift_soln(c, dtuse, L, 2, X_edge, Z_edge, uedge, Z_newgrid_edge);

if(iso==0 || layers_add_bottom==0)
    wcent_t1 = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, wcent, Z_st);
    wcent_t2 = timeshift_soln(c, dtuse, L, 2, XC_st, Z_st, wcent, Z_st);
else
wcent_t1 = timeshift_soln(c, dtuse, L, 1, X_shift_w, Z_shift_w, wcent, Z_shift_w);
wcent_t2 = timeshift_soln(c, dtuse, L, 2, X_shift_w, Z_shift_w, wcent, Z_shift_w);
end

% wcent_t1 = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, wcent, Z_newgrid);
% wcent_t2 = timeshift_soln(c, dtuse, L, 2, XC_st, Z_st, wcent, Z_newgrid);
% wedge_t1 = timeshift_soln(c, dtuse, L, 1, X_edge, Z_edge, wedge);
% wedge_t2 = timeshift_soln(c, dtuse, L, 2, X_edge, Z_edge, wedge);

density_t1 = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, density, Z_newgrid);
density_t2 = timeshift_soln(c, dtuse, L, 2, XC_st, Z_st, density, Z_newgrid);

%q_full = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, q_layer_av, Z_newgrid);
q_full = pnhz;

if(iso)
    dzz_t1 = timeshift_soln(c, dtuse, L, 1, XC_st, Z_st, layer_thickness, Z_newgrid);
    dzz_t2 = timeshift_soln(c, dtuse, L, 2, XC_st, Z_st, layer_thickness, Z_newgrid);
    dzz_t3 = timeshift_soln(c, dtuse, L, 3, XC_st, Z_st, layer_thickness, Z_newgrid);
    q_full = q_layer_av; 
end

%%


if(iso)
 Z_w_st = cumsum(layer_thickness, 1);
 Z_w_st = cat(1, zeros(1, length(xc_st)), Z_w_st);
 Z_w_st = Z_w_st + -300;

  Z_w_st_t1 = cumsum(dzz_t1, 1);                                                                                                                                
 Z_w_st_t1 = cat(1, zeros(1, length(xc_st)), Z_w_st_t1);                                                                                                                   
 Z_w_st_t1 = Z_w_st_t1 + -300;

 Z_w_st_t2 = cumsum(dzz_t2, 1);                                                                                                                                
 Z_w_st_t2 = cat(1, zeros(1, length(xc_st)), Z_w_st_t2);                                                                                                                   
 Z_w_st_t2 = Z_w_st_t2 + -300;

% sherlock just had these three lines commneted in, rest commented out
 % Z_other_st = cumsum(layer_thickness, 1) - layer_thickness./2 - 300; 
 % Z_other_st_t1 = cumsum(dzz_t1, 1) - dzz_t1./2 - 300;                                                                                                 
 % Z_other_st_t2 = cumsum(dzz_t2, 1) - dzz_t2./2 - 300;  

 Z_other_st(1, :) = -300+layer_thickness(1, :)./2; 
 Z_other_st_t1(1, :) = -300+dzz_t1(1, :)./2;                                                                                  
 Z_other_st_t2(1, :) = -300+dzz_t2(1, :)./2; 
 for k = 2:nk_cont
     Z_other_st(k, :) = Z_other_st(k-1, :)+layer_thickness(k-1, :)./2+layer_thickness(k, :)./2; 
    Z_other_st_t1(k, :) = Z_other_st(k-1, :)+dzz_t1(k-1, :)./2+dzz_t1(k, :)./2;                                                                                  
    Z_other_st_t2(k, :) = Z_other_st(k-1, :)+dzz_t1(k-1, :)./2+dzz_t2(k, :)./2; 
 end

 Z_hyrbid_w = cat(1, Z_w_st(1:layers_add_bottom+1, :),  Z_other_st(layers_add_bottom+1:end, :));
 Z_hyrbid_w_t1 = cat(1, Z_w_st_t1(1:layers_add_bottom+1, :),  Z_other_st_t1(layers_add_bottom+1:end, :));
 Z_hyrbid_w_t2 = cat(1, Z_w_st_t2(1:layers_add_bottom+1, :),  Z_other_st_t2(layers_add_bottom+1:end, :));


else
 Z_w_st = cumsum(dz_st(end:-1:1));
 Z_w_st = repmat(Z_w_st', length(xv), 1);
 Z_w_st = Z_w_st';

 Z_w_st = cat(1, zeros(1, length(xv)), Z_w_st);
 Z_w_st = Z_w_st + -300;
 Z_w_st_t1 = Z_w_st;
 Z_w_st_t2 = Z_w_st;

 Z_other_st = cumsum(dz_st(end:-1:1)) - dz_st(end:-1:1)./2 - 300;
 Z_other_st = repmat(Z_other_st', length(xv), 1);
 Z_other_st = Z_other_st';
 % [~, Z_other_st] = meshgrid(xc_st, ze);
 nk_cont = Nk;
 Z_hyrbid_w = Z_st;
end

wcent_new = zeros(nk_cont+1, length(xc_st));
wcent_newt1 = zeros(nk_cont+1, length(xc_st));
wcent_newt2 = zeros(nk_cont+1, length(xc_st));

XC_wst = cat(1, XC_st(1, :),  XC_st);

if(iso ==0 || layers_add_bottom==0)
    Z_hyrbid_w = Z_other_st;
    XC_wst_old = XC_st;
else
    XC_wst_old = XC_wst;
end

for i=1:length(xc_st) 
    wcent_new(2:end-1, i) = interp1(Z_hyrbid_w(:, i), wcent(:, i), Z_w_st(2:end-1, i));
    wcent_newt1(2:end-1, i) = interp1(Z_hyrbid_w(:, i), wcent_t1(:, i), Z_w_st_t1(2:end-1, i));
    wcent_newt2(2:end-1, i) = interp1(Z_hyrbid_w(:, i), wcent_t2(:, i), Z_w_st_t2(2:end-1, i));

end
%XC_wst = cat(1, XC_st(1, :),  XC_st);

wcent_old = wcent;
%figure;
%pcolor(XC_wst, Z_hyrbid_w, wcent_old);
%colormap jet;
%shading interp;


%figure;
%pcolor(XC_wst, Z_w_st, wcent_new);
%colormap jet;
%shading interp;


%% 
wcent = wcent_new;
wcent_t1 = wcent_newt1;
wcent_t2 = wcent_newt2;



%% calculate richardson number? 

%say dudx is u(up)-u(down)/dz
% 
% dudz = zeros(nk_cont-1, Nx);
% Z_ri = zeros(nk_cont+1, Nx);
% Z_ri(1, :) = -300;
% for k = 1:nk_cont-1
%     dudz(k, :) = (ucent(k, :)-ucent(k+1, :))/(.5.*(layer_thickness(k, :)+layer_thickness(k+1, :)));
%     Z_ri(k+1, :) = Z_ri(k, :) + layer_thickness(k+1, :);
% end
% dudz = dudz.^2;
% 
% Ri = N2(Z_ri(2:end-1, :))./dudz;
% 
% figure;
% pcolor(XC_st(2:end, :), Z_ri(2:end-1, :), Ri)
% shading flat
% colorbar
% 


%%
%dt = 20; %set arbitrary dt for testing

%Z_st = cumsum(layer_thickness, 1);

% % %shift everything back one and two time steps to write t0-dt and t0-2dt
% ucent_t1 = timeshift_soln(c, dt, L, 1, XC_st, Z_st, ucent);
% ucent_t2 = timeshift_soln(c, dt, L, 2, XC_st, Z_st, ucent);
% uedge_t1 = timeshift_soln(c, dt, L, 1, X_edge, Z_edge, uedge);
% uedge_t2 = timeshift_soln(c, dt, L, 2, X_edge, Z_edge, uedge);
% 
% wcent_t1 = timeshift_soln(c, dt, L, 1, XC_st, Z_st, wcent);
% wcent_t2 = timeshift_soln(c, dt, L, 2, XC_st, Z_st, wcent);
% % wedge_t1 = timeshift_soln(c, dt, L, 1, X_edge, Z_edge, wedge);
% % wedge_t2 = timeshift_soln(c, dt, L, 2, X_edge, Z_edge, wedge);
% 
% density_t1 = timeshift_soln(c, dt, L, 1, XC_st, Z_st, density);
% density_t2 = timeshift_soln(c, dt, L, 2, XC_st, Z_st, density);
% 
% if(iso)
%     dzz_t1 = timeshift_soln(c, dt, L, 1, XC_st, Z_st, layer_thickness);
%     dzz_t2 = timeshift_soln(c, dt, L, 2, XC_st, Z_st, layer_thickness);
% end

%Z_st = cumsum(layer_thickness, 1);

if(iso)
    
    min_dzz = min(min(layer_thickness))
    min_dzz_hybrid =  min(min(layer_thickness(Nkbot+1:end, :)))
    
    layers_add_mid
end

% 
%shading flat
% 
% figure(3);
% pcolor(XC_st, Z_st, density_t2);
% %shading flat
% 
% figure(4)
% pcolor(XC_st, Z_st, abs(density-density_t1));
% shading flat
% colorbar
% 
% figure(5)
% pcolor(XC_st, Z_st, abs(ucent-ucent_t2));
% shading flat
% colorbar

%%
close all;


dt = getvalue(suntansfile,'dt');
dtuse = dt;

Nsteps = L/c * 10 / dt 

Courantnum = dt*c/dx
Courantnum_alt = dt*max(max(abs(ucent)))/dx %why would it be u+c?


%%
%%% Output for input into SUNTANS %%%
numprocs=0;
while(exist([datadir,'/edgedata.dat.',num2str(numprocs)],'file'))
    numprocs=numprocs+1;
end
fprintf('Found %d processors.\n',numprocs);

if(WRITE)
    for proc=0:numprocs-1

        % Load the edge-centered grid data
        edges = load([datadir,'/edgedata.dat.',num2str(proc)]);

        % Components of vector normal to edges
        n1 = edges(:,3)*ones(1,Nk);
        n2 = edges(:,4)*ones(1,Nk);
        % Centers of edges
        xe = edges(:,5);

        % Number of  edges
        Ne = length(xe);

        % Create the X,Z grid on edge centers
        [Xe,Ze]=ndgrid(xe,1:Nk);

        %match of velocities with grid for each processor
        Ue = zeros(Ne, Nk);
        Ue_t1 = zeros(Ne, Nk);
        Ue_t2 = zeros(Ne, Nk);

        tol = dx/100; 

        for xval = 1:Ne
	    ind = find(abs(xe_st - xe(xval))<tol);
            if isempty(ind)
                centind = find(abs(xc_st - xe(xval))<tol);
                Ue(xval, :) = ucent(:, centind);
                Ue_t1(xval, :) = ucent_t1(:, centind);
                Ue_t2(xval, :) = ucent_t2(:, centind);

            else
                Ue(xval, :) = uedge(:, ind);
                Ue_t1(xval, :) = uedge_t1(:, ind);
                Ue_t2(xval, :) = uedge_t2(:, ind);

                % We_t1(xval, :) = wedge_t1(:, ind);
                % We_t2(xval, :) = wedge_t2(:, ind);
                % We(xval, :) = wedge(:, ind);
            end 

        end

        Ve = zeros(Ne,Nk);
        if nudge==1
           Ue = Ue - c; %subtract off c to stay inline with wave
	    end
        %Ue_t1 = Ue_t1 - c; %subtract off c to stay inline with wave
        %Ue_t2 = Ue_t2 - c; %subtract off c to stay inline with wave


        % U is the dot product of the internal wave velocity with the normal
        U = Ue.*n1 + Ve.*n2;
        U_t1 = Ue_t1.*n1 + Ve.*n2;
        U_t2 = Ue_t2.*n1 + Ve.*n2;

        % Load the cell-centered data
        cells = load([datadir,'/celldata.dat.',num2str(proc)]);

        % Cell centers
        xv = cells(:,2);

        % Number of cells
        Nc = length(xv);

        %try guessing h 
        a1 = 1.765e-6;
        b1 = 5000;
        c1 = 536.6;
        min_h = -1.9417e-7;

        hfunc = @(xval) min_h + a1*exp(-((xval-b1)/c1)^2); 
        hguess = zeros(Nc, 1);

        % W_t1 = zeros(Nc, Nk);
        % W_t2 = zeros(Nc, Nk);
        % W = zeros(Nc, Nk);

        % Create the X,Z grid on cell centers
        [Xv,Zv]=ndgrid(xv,1:Nk);

        % Salinity and temperature are defined on cell centers
        S_full = density - rho0; %density perturbation 
        S_full_t1 = density_t1 - rho0; %density perturbation 
        S_full_t2 = density_t2 - rho0; %density perturbation 

        indskeep = [];
        for xval = 1:Nc
            %ind = find(xc == xv(xval));
            ind = find(abs(xc_st - xv(xval))<tol);
            indskeep = [indskeep, ind];
            hguess(xval) = hfunc(xv(xval));
        end

        W_t1 = wcent_t1(:, indskeep);
        W_t2 = wcent_t2(:, indskeep);
        W = wcent(:, indskeep);

        W_t1=W_t1';
        W_t2=W_t2';
        W=W';

        S = S_full(:, indskeep);
        S = S';

        S_t1 = S_full_t1(:, indskeep);
        S_t1 = S_t1';

        S_t2 = S_full_t2(:, indskeep);
        S_t2 = S_t2';
        if(iso)
            layer_thickness_new = layer_thickness(:, indskeep);
            layer_thickness_new = layer_thickness_new';

            layer_thickness_t1 = dzz_t1(:, indskeep);
            layer_thickness_t1 = layer_thickness_t1';

            layer_thickness_t2 = dzz_t2(:, indskeep);
            layer_thickness_t2 = layer_thickness_t2';

            layer_thickness_t3 = dzz_t3(:, indskeep);
            layer_thickness_t3 = layer_thickness_t3';
        end
        %S = 0*S; %set s = 1 as a test

	    q = q_full(:, indskeep);
        q=q';
        fs = eta(end, indskeep); %use as free surface initialization

        %put T=1 in non isopyncals, 0 in isopycnals
        T = zeros(Nc,Nk);
        % if(iso)
        %     %T(:, 1:layers_add_bottom+layers_add_mid+layers_add_top)=1;
        %     T(xv>4500 & xv<5500, :) = 1;
        % else
        %     T(xv>4500 & xv<5500, :) = 1;
        %     %T = T + 1;
        % end
        %T = ones(Nc,Nk);


        %T(:, 1:3) = 1; %put T=1 in bottom three layers

        % Initial free surface is zero (not currently in use)
        h = zeros(Nc,1);

        u_file = [u_init_file,'.',num2str(proc)];
        w_file = [w_init_file,'.',num2str(proc)];
        s_file = [s_init_file,'.',num2str(proc)];
        T_file = [T_init_file,'.',num2str(proc)];
        fs_file = [fs_init_file,'.',num2str(proc)];   
        dzz_file = [dzz_init_file,'.',num2str(proc)];
        dzz_t1_file = [dzz_t1_init_file,'.',num2str(proc)];
        dzz_t2_file = [dzz_t2_init_file,'.',num2str(proc)];
        dzz_t3_file = [dzz_t3_init_file,'.',num2str(proc)];
        u_fid = fopen(u_file,'wb');
        w_fid = fopen(w_file,'wb');
        s_fid = fopen(s_file,'wb');
        T_fid = fopen(T_file,'wb');
        fs_fid = fopen(fs_file,'wb');
        dzz_fid = fopen(dzz_file,'wb');
        dzz_t1_fid = fopen(dzz_t1_file,'wb');
        dzz_t2_fid = fopen(dzz_t2_file,'wb');
        dzz_t3_fid = fopen(dzz_t3_file,'wb');

        %previous time steps
        u_t1_fileid = [u_t1_file,'.',num2str(proc)];
        u_t2_fileid = [u_t2_file,'.',num2str(proc)];
        w_t1_fileid = [w_t1_file,'.',num2str(proc)];
        w_t2_fileid = [w_t2_file,'.',num2str(proc)];
        s_t1_fileid = [s_t1_file,'.',num2str(proc)];    
        s_t2_fileid = [s_t2_file,'.',num2str(proc)];   
        q_fileid = [q_init_file,'.',num2str(proc)];
        u_t1_fid = fopen(u_t1_fileid,'wb');
        u_t2_fid = fopen(u_t2_fileid,'wb');
        w_t1_fid = fopen(w_t1_fileid,'wb');
        w_t2_fid = fopen(w_t2_fileid,'wb');
        s_t1_fid = fopen(s_t1_fileid,'wb');
        s_t2_fid = fopen(s_t2_fileid,'wb');
        q_fid = fopen(q_fileid,'wb');

        % Write the data, noting that suntans expects z values
        % decreasing with increasing index

        mult0 = 1;

        for j=1:Ne
            fwrite(u_fid,U(j,end:-1:1)*mult0,'float64');
            
            fwrite(u_t1_fid,U_t1(j,end:-1:1)*mult0,'float64');
            fwrite(u_t2_fid,U_t2(j,end:-1:1)*mult0,'float64');

        end
        for i=1:Nc
            fwrite(s_fid,S(i,end:-1:1)/rho0*mult0,'float64');
            fwrite(s_t1_fid,S_t1(i,end:-1:1)/rho0*mult0,'float64');
            fwrite(s_t2_fid,S_t2(i,end:-1:1)/rho0*mult0,'float64');
            fwrite(T_fid,T(i,end:-1:1),'float64'); 
            fwrite(q_fid,q(i,end:-1:1),'float64'); 
            fwrite(w_fid,W(i,end:-1:1)*mult0,'float64');
            fwrite(w_t1_fid,W_t1(i,end:-1:1)*mult0,'float64');
            fwrite(w_t2_fid,W_t2(i,end:-1:1)*mult0,'float64');
            if(iso)
                fwrite(dzz_fid,layer_thickness_new(i,end:-1:1),'float64');
                fwrite(dzz_t1_fid,layer_thickness_t1(i,end:-1:1),'float64');
                fwrite(dzz_t2_fid,layer_thickness_t2(i,end:-1:1),'float64');
                fwrite(dzz_t3_fid,layer_thickness_t3(i,end:-1:1),'float64');
            end
            
        end


        fwrite(fs_fid,h,'float64'); %if rigid lid
        %fwrite(fs_fid,fs','float64');
        %fwrite(fs_fid, hguess, 'float64');
        

        %fprintf('Created %s, %s, %s, %s, %s, and %s.\n',u_file,u_t1_fileid,u_t2_fileid,s_file,T_file,fs_file);
        fclose(u_fid);
        fclose(w_fid);
        fclose(s_fid);
        fclose(T_fid);
        fclose(dzz_fid);
        fclose(dzz_t1_fid);
        fclose(dzz_t2_fid);
        fclose(dzz_t3_fid);
        fclose(fs_fid);

        fclose(u_t1_fid);
        fclose(u_t2_fid);
        fclose(w_t1_fid);
        fclose(w_t2_fid);
        fclose(s_t1_fid);
        fclose(s_t2_fid);
        fclose(q_fid);

    end
end
