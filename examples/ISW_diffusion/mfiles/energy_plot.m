
clear all; close all;


ANIMATE = ~false;



datadir='../data_diffusion_test';

addpath '/Users/brookepauken/suntans/suntans-master/mfiles';


suntansfile=[datadir,'/suntans.dat'];
 
nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
beta = getvalue(suntansfile,'beta');
dt = getvalue(suntansfile,'dt');
nout = 1+nsteps/ntout;

cwave = 1.76776;

g = 9.81;
rho0 = 1000;

[x,z,s]=plotslice('s',datadir,1);
[Nz,Nx]=size(x);
L = 10000;
H = 300;
dx = L/Nx;
dz = H/Nz;

Ep = zeros(nout,1);
Ek = zeros(nout,1);
Eb = zeros(nout,1);

norm_rho = zeros(nout,1);
norm_rhodv = zeros(nout,1);
norm_N2 = zeros(nout,1);
norm_eta = zeros(nout,1);
max_eta = zeros(nout,1);

depthdata = load([datadir,'/depth.dat-voro'],'-ascii');
dv = depthdata(:,3);

iso = 1;

targetgrid='interior';

figure(1);
clf;
for n=1:1:nout   

    [x,z,s]=plotslice('s',datadir,n);
    [x,z,u]=plotslice('u',datadir,n);
    [x,z,w]=plotslice('w',datadir,n);

    if(n==1)
        s0 = s;
        x0 = x;
        z0 = z;
        s0 = rho0+s0*rho0;

        Nkmax = getvalue(suntansfile,'Nkmax');
        Nc = length(dv);

        fid = fopen([datadir,'/vert_dzz.dat'],'rb');
        fseek(fid,8*Nkmax*Nc*(n-1),'bof');
        dzz = reshape(fread(fid,Nkmax*Nc,'float64'),Nc,Nkmax);
        dzz = dzz';
        fclose(fid);

        s0 = s0;
        u0 = sqrt(u.^2 + w.^2);
        dzz0 = dzz;


        [drho_dz, N20] = computeN2(s0, z, Nkmax, Nc);
        
        snew = s0;
        dzznew = dzz0;
        N2new = N20;

        rhodv0 = (snew.*dzznew.*dx);

        z_bar = linspace(0, -300, 1000);
        [rho_background, rho_solve, rho] = get_backgroundrho(z_bar);

        z_solve = z;

        initetaT = tic;

        for i=1:Nc
            for k=1:Nkmax
                specific_rhosolve = @(z) rho_solve(z,s0(k, i));
                %specific_rhosolve(z(i,k))
                z_solve(k, i) = fzero(specific_rhosolve, 0);
            end
        end

        Tprint = toc(initetaT);
        fprintf('Inital get eta took %8.2f seconds \n', Tprint);


        eta0 = z-z_solve;
        etanew = eta0;

    else

        %at n = 1, timestep is dt
        %at n = 2, timestep is dt + dt*ntout
        %at n = 3, timestep is dt + dt*ntout*2

        %so want to shift it ntout*(n-1)
        zshift = 1:Nkmax;
        zshift = repmat(zshift,1,Nc);

        if(iso)
            if(n==2)
               snew = timeshift_soln(cwave, dt, 10000, 1, x, zshift, s0, x);
               dzznew = timeshift_soln(cwave, dt, 10000, (1), x, zshift, dzz0, x);
            else
                snew = timeshift_soln(cwave, dt, 10000, (n)*ntout, x, zshift, s0, x);
                dzznew = timeshift_soln(cwave, dt, 10000, (n)*ntout, x, zshift, dzz0, x);
            end


            %if n==3, want to shift 1
            snew = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, s0, x);
            dzznew = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, dzz0, x);
            N2new = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, N20, x);
            etanew = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, eta0, x);
    
        else
            snew = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, s0, x);
            dzznew=dzz0;
            N2new = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, N20, x);
            etanew = timeshift_soln(cwave, dt, 10000, (n-1)*ntout, x, zshift, eta0, x);
        end


    end


    Nkmax = getvalue(suntansfile,'Nkmax');
    Nc = length(dv);

    NX = Nc;
    NZ = Nkmax;
    Ubg = @(z) z*0;
    Ubgz = @(z) z*0;

    fid = fopen([datadir,'/vert_dzz.dat'],'rb');
    fseek(fid,8*Nkmax*Nc*(n-1),'bof');
    dzz = reshape(fread(fid,Nkmax*Nc,'float64'),Nc,Nkmax);
    dzz = dzz';
    fclose(fid);
    
    s = rho0+s*rho0;
    [s_sorted, inds] = sort(s(:));
    dz_sorted = dx*dzz(inds)/L;

    [drho_dz, N2] = computeN2(s, z, Nkmax, Nc);


    z_solve = z;
    for i=1:Nc
        for k=1:Nkmax
            specific_rhosolve = @(z) rho_solve(z,s(k, i));
            %specific_rhosolve(z(i,k))
            z_solve(k, i) = fzero(specific_rhosolve, z(k,i)-etanew(k,i));
        end
    end



    eta = z-z_solve;

    z_sorted = -cumsum(dz_sorted);
    z_sorted = [0; z_sorted];
    z_sorted = .5*(z_sorted(2:end, :)+z_sorted(1:end-1, :));


    Eb(n) = beta*sum(s_sorted(:).*z_sorted(:)*g.*dz_sorted(:))*L;
    Ep(n) = beta*sum(sum(s.*z.*g.*dzz))*dx;
    Ek(n) = 0.5*sum(sum(s.*(u.^2+w.^2).*dzz))*dx;
    norm_rho(n) = norm(s-snew);
    norm_rhodv(n) = norm(((s-rho0).*dzz.*dx-(snew-rho0).*dzznew.*dx));
    norm_N2(n) = norm(N2-N2new);
    norm_eta(n) = norm(eta-etanew);
    % norm_eta(n) = norm(residual);
    max_eta(n) = max(max(abs(eta)))-max(max(abs(eta0)));

    if(ANIMATE)
        tiledlayout("vertical");

        zplot_current(1, :) = -dzz(1, :)/2;
        zplot_new(1, :) = -dzznew(1, :)/2;
        for(k=2:Nkmax)
            zplot_current(k, :) = zplot_current(k-1, :) - (dzz(k-1, :)./2 + dzz(k, :)./2);
            zplot_new(k, :) = zplot_new(k-1, :) - (dzznew(k-1, :)./2 + dzznew(k, :)./2);
        end

        %xplot_interp = meshgrid(x(1, :), z_cartesian(:, 1));

        nexttile
        pcolor(x,zplot_current,eta);
        colorbar
        shading flat

        nexttile
        pcolor(x,zplot_new,etanew);
        colorbar
        shading flat


        nexttile;
        pcolor(x,zplot_new,(eta-etanew));
        colorbar;
        %cmocean('balance', 'pivot',0)
        shading flat;

        %axis image;
        nexttile;
        pcolor(x,zplot_new,dzz-dzznew);
        colorbar;
        %cmocean('balance', 'pivot',0)
        shading flat;
        %axis image;
        
        drawnow;

    end
end

%%
figure;
t = [1:nout]*ntout*dt;
plot(t, norm_rho)

figure;
plot(t, norm_rhodv)

figure;
plot(t, norm_N2)



if(n~=nout)
    nout = n-1;
end

Eb = Eb(1:nout);
Ep = Ep(1:nout);
Ek = Ek(1:nout);
Ea = Ep-Eb;

Ep0 = Ep(1);
Et = Ek+Ep-Ep0;
t = [1:nout]*ntout*dt;

etamix = (Eb(1)-Eb(end))/Eb(1);



pline = '-';
fprintf('Mixing efficiency = %.2f\n',etamix);
Enorm = abs(Ep0);
figure(2)
plot(t,Eb/Enorm,['r',pline],t,Ek/Enorm,['b',pline],t,(Ep-Ep0)/Enorm,['g',pline],t,Et/Enorm,['k',pline]);
xlabel('t (s)');
ylabel('Energy');
legend('E_b/E_{p0}','E_k/E_{p0}','(E_p-E_{p0})/E_{p0}','E_t/E_{p0}');
set(gca,'fontsize',14);
grid on

figure;
%plot(t,Ek-Ek(1))
hold on
plot(t, Ep-Ep(1))
hold on
%plot(t, Eb-Eb(1))
plot(t, Ea-Ea(1))
plot(t, Et-Et(1))
legend('Ep', 'Ea', 'Et')

figure;
plot(t,Ek-Ek(1))
hold on
plot(t,Ea-Ea(1))
hold on
plot(t,Et-Et(1))
plot(t, Eb-Eb(1))
legend('Ek', 'Ea', 'Et', 'Eb')



%% FUNCTION 

%Compute N2

function [drho_dz, N2] = computeN2(rho_bar, z, Nk, Nc)
    g = 9.81;
    rho0 = 1000;

    drho_dz = zeros(Nk,Nc);
    for k=2:Nk-1
        drho_dz(k, :) = (rho_bar(k+1, :)-rho_bar(k-1, :))./(z(k+1, :)-z(k-1, :));
    end
    drho_dz(1, :) = 1.5*drho_dz(2, :)-0.5*drho_dz(3, :);
    drho_dz(Nk, :) = 1.5*drho_dz(Nk-1, :)-0.5*drho_dz(Nk-2, :);
    N2 = (-g/rho0*drho_dz);
end

%Compute background density at more dense z 
function [rho_background, rho_solve, rho] = get_backgroundrho(z_bar)
    % Stratification is given by a tanh distribution
    D = 300;         % Depth
    multiplier = 0.25;
    h1 = multiplier*D; %0.005*D; %0.25*D;  (.1)       % Height of pycnocline
    Delta_rho = 6; %6;       % Density difference
    rho0 = 1000;         % Reference density
    g = 9.81;
    delta = 0.15*D; %0.1*D;       % Pycnocline thickness
    
    lin_start = -180; %was -180;
    lin_start = lin_start/300*D;


    if(1)
        %lin_start = -301; %to turn off
        rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))); 
        rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2);
       
       
        %continue slope below lin_start down
        m_lin = rhoz(lin_start);
    end


    if(0)
        rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))+rho0); 
        rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2);
    else
        rho = @(z) (Delta_rho/2*(1-tanh((z+h1)/delta))+rho0).*(z>=lin_start) + (m_lin*(z-lin_start) + rho(lin_start) + rho0).* (z<lin_start); 
        rhoz= @(z) -((Delta_rho/(2*delta))*sech((z+h1)/delta).^2).*(z>=lin_start) + (m_lin).* (z<lin_start);
        rho_solve = @(z, s_input) rho(z) - s_input;
    end

    rho_background = rho(z_bar);
end


function [residual, LHS, RHS] = calc_residual(ks, ms, eta, Ubg, Ubgz, N2, z, c, gridtype)
    % DJL residual using Eq 2.32 in (Stastna, 2001)
    % Adapted from DJLES code for use here
    
    % Odd extend the function in both directions
    etaextended = djles_extend(eta, 'odd', 'odd', gridtype);
    
    [SZ,SX] = size(eta);
    
    % Compute left hand side
    LAP = -bsxfun(@plus,ms.^2,ks.^2);           % Laplacian operator
    LHS = real(ifft2(LAP.*fft2(etaextended)));  % Laplacian of extended eta
    LHS = LHS(1:SZ, 1:SX);                      % Trim for eta on gridtype
    
    % Compute right hand side
    [etax, etaz] = djles_gradient(eta, ks, ms, 'odd', 'odd', gridtype);
    Umc = Ubg(z-eta)-c;
    RHS = -(Ubgz(z-eta)./Umc).*(1 - (etax.^2 + (1-etaz).^2)) - N2.*eta./(Umc.^2);
    
    % Residual
    residual = LHS-RHS;
end

